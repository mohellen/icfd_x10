package nssim;
import x10.io.File;
import x10.io.FileWriter;
import x10.io.FileReader;
import x10.io.Printer;
import x10.compiler.Inline;
import x10.util.List;
import x10.util.ArrayList;

public class NSSim extends Parameters{

  /**********************
  *  Constructor
  * *********************/
  public def this()
  {
    super();
  _ncx = _ncx * _resx;
  _ncy = _ncy * _resy;
  _dx = _domain_size_x / _ncx;
  _dy = _domain_size_y / _ncy;

  }//end constructor

  public def run(): void
  {
    var dt: Double = 0.0;
    var t:  Double = 0.0;
    val t_end = final_time;
    var time_step: Int = 0;

    var vtk_freq:Double = 0.1;
    var vtk_cnt:Int = 1;
    val vtk_outfile:String = output_name;

    // Create global mask and global fields
    val global_mask = GlobalRef[Array[Int]{rank==2}](create_geometry_mask());
    val global_u = GlobalRef[Array[Double]{rank==2}](new Array[Double](0..((_ncy+1)) * 0..((_ncx+1)), _initial_velocity_x));
    val global_v = GlobalRef[Array[Double]{rank==2}](new Array[Double](0..((_ncy+1)) * 0..((_ncx+1)), _initial_velocity_y));
    val global_p = GlobalRef[Array[Double]{rank==2}](new Array[Double](0..((_ncy+1)) * 0..((_ncx+1)), _initial_pressure));
    val global_f = GlobalRef[Array[Double]{rank==2}](new Array[Double](0..((_ncy+1)) * 0..((_ncx+1)), 0.0));
    val global_g = GlobalRef[Array[Double]{rank==2}](new Array[Double](0..((_ncy+1)) * 0..((_ncx+1)), 0.0));

    val map: Array[Int] = new Array[Int](1 , 0);
    var b:NSBlock = create_block(map , 0, 1, _ncx, _ncy, global_mask, global_u, global_v, global_p, global_f, global_g);

    // Initialize U,V,P,F,G
    b.update_domain_uv();
    b.update_domain_p();
    b.update_domain_fg();

    // Initial Output
    if (vtk_flag==1) {
      b.write_vts_file(vtk_outfile, map, 0, 1);
      b.write_vtm_file(vtk_outfile, 0, 1);
    }


    // Main loop
    while (t < t_end) {
      // 1. compute dt
      b.compute_local_dt();

      // 2. compoute F,G
      b.compute_fg();
      b.update_domain_fg();

      // 3. solve for pressure P
      b.compute_p();
      b.update_domain_p();

      // 4. Compute velocity U,V
      b.compute_uv();
      b.update_domain_uv();
      b.update_boundary_uv();

      // 5. Update simulation time
      t += b.get_dt();
      Console.OUT.printf("Simulation completed at t = %.5f, dt = %.5f \n", t, b.get_dt());

      // 6. Write VTK output
      if (vtk_flag==1) {
        if ((Math.floor(t/vtk_freq) as Int) >= vtk_cnt) {
          Console.OUT.printf("Writing output files at t = %.5f\n", t);
          b.write_vts_file(vtk_outfile, map, vtk_cnt, 1);
          b.write_vtm_file(vtk_outfile, vtk_cnt, 1);
          vtk_cnt=vtk_cnt+1;
        }
      }

      time_step+=1;
      Console.OUT.printf("at the end of time step %d, num proc = %d \n", time_step, 1);
    }//end while

    return;
  }// end funcion run()

  /********************
  *  Internal Methods
  ********************/
  protected def create_block(
  map: Array[Int]{rank==1},
  var blockid: Int,  // The block ID
  nblocks: Int,  // Total number of blocks (desired)
  gncx: Int, 	   // Global resolution number of inner cells in x-dimension
  gncy: Int,     // Global resolution number of inner cells in x-dimension
  global_mask: GlobalRef[Array[Int]{rank==2}],
  global_u: GlobalRef[Array[Double]{rank==2}],
  global_v: GlobalRef[Array[Double]{rank==2}],
  global_p: GlobalRef[Array[Double]{rank==2}],
  global_f: GlobalRef[Array[Double]{rank==2}],
  global_g: GlobalRef[Array[Double]{rank==2}]
  ) : NSBlock
  {

    //map inverse
    for (i in (0..(nblocks-1))  ){
      if (map(i)==blockid)
      blockid = i;
    }

    // We only use to 2^n blocks
    val n = Math.floor(Math.log(nblocks) / Math.log(2)) as Int;
    val nbx: Int;  // number of blocks in x-direction
    val nby: Int;  // number of blocks in y-direction

    // Grid of Blocks: nbx*nby, row-major
    if ((n%2) as Int == 0) {
      nbx = Math.pow(2,n/2) as Int;
      nby = Math.pow(2,n/2) as Int;
    } else {
      nbx = Math.pow(2,n/2+1) as Int;
      nby = Math.pow(2,n/2) as Int;
    }
    val nb = nbx*nby; // the actual number of blocks

    // Warning if actual # blocks is smaller than input
    if (nb < nblocks) {
      Console.OUT.printf("The actual number of blocks %d is smaller than the desired number %d.\n", nb, nblocks);
    }

    // Compute neighbor block IDs
    // A negative neighbor id indicates domain boundary (no neighbor)
    var nid: Int = ((blockid + nbx) < nb) ? (blockid + nbx) : (-1);
    var sid: Int = ((blockid - nbx) >= 0) ? (blockid - nbx) : (-1);
    var eid: Int = ((blockid + 1)% nbx == 0) ? (-1) : (blockid + 1);
    var wid: Int = (blockid % nbx == 0) ? (-1) : (blockid - 1);

    // Compute block index range
    var gimin: Int;	// Glboal index of local inner domain imin (x-dimension, inner loop)
    var gimax: Int; // Glboal index of local inner domain imax
    var gjmin: Int; // Glboal index of local inner domain jmin (y-dimension, outer loop)
    var gjmax: Int; // Glboal index of local inner domain jmax

    val seq = (bj:Int, bi:Int) => {bj*nbx + bi}; // conversion of 2D index to 1D index
    val bi = (blockid % nbx) as Int;
    val bj = (blockid / nbx) as Int;
    val xtrunk = (gncx / nbx) as Int;
    val ytrunk = (gncy / nby) as Int;

    gimin = xtrunk*bi + 1;
    gimax = gimin + xtrunk - 1;
    if ((bi+1)%nbx == 0) gimax += gncx % nbx;

    gjmin = ytrunk*bj + 1;
    gjmax = gjmin + ytrunk - 1;
    if ((bj+1)%nby == 0) gjmax += gncy % nby;

    // map the ids to the new id
    blockid=map(blockid);

    if (nid !=-1)
    nid=map(nid);
    if (sid !=-1)
    sid=map(sid);
    if (eid !=-1)
    eid=map(eid);
    if (wid !=-1)
    wid=map(wid);

    return new NSBlock( blockid, nid, sid, eid, wid,
    gimin, gimax, gjmin, gjmax, gncx, gncy, global_mask, global_u, global_v, global_p, global_f, global_g);
  }//end function

  protected def create_geometry_mask(): Array[Int]{rank==2}
  {
    /************************************************
    *  Fluid cells are marked with "0".
    *  Non-fluid cells are marked with "1xxxx".
    *
    *  Each non-fluid cells contains 5-bit:
    *      [ C  E  W  S  N ]
    *   C: center, 		cell [j  ][i  ]
    *   E: east  neighbor, cell [j  ][i+1]
    *   W: west  neighbor, cell [j  ][i-1]
    *   N: north neighbor, cell [j+1][i  ]
    *   S: south neighbor, cell [j-1][i  ]
    *
    *  "1" indicates an obstacle cell
    *  "0" indicates a fluid cell
    * **********************************************/
    val mask = new Array[Int](0..((_ncy+1)) * 0..((_ncx+1)), FLUID);

    // 1. Initialize all non-fluid cells to 10000
    // 1.1 Setup boundaries
    if (_boundary_west == BOUNDARY_TYPE_NOSLIP || _boundary_west == BOUNDARY_TYPE_FREESLIP) {
      for (j in 1.._ncy)
      mask(j,0) = 10000;
    }
    if (_boundary_east == BOUNDARY_TYPE_NOSLIP || _boundary_east == BOUNDARY_TYPE_FREESLIP) {
      for (j in 1.._ncy)
      mask(j,_ncx+1) = 10000;
    }
    if (_boundary_north == BOUNDARY_TYPE_NOSLIP || _boundary_north == BOUNDARY_TYPE_FREESLIP) {
      for (i in 0..(_ncx+1))
      mask(_ncy+1,i) = 10000;
    }
    if (_boundary_south == BOUNDARY_TYPE_NOSLIP || _boundary_south == BOUNDARY_TYPE_FREESLIP) {
      for (i in 0..(_ncx+1))
      mask(0,i) = 10000;
    }

    // 1.2 Setup inner domain
    var ncl: Int;
    var nch: Int;
    var nrl: Int;
    var nrh: Int;
    for (k in 0..(NUM_OBS-1)) {
      ncl = (Math.round(OBS_XMIN(k)/_dx) as Int) + 1;
      nch = (Math.round(OBS_XMAX(k)/_dx) as Int);
      nrl = (Math.round(OBS_YMIN(k)/_dy) as Int) + 1;
      nrh = (Math.round(OBS_YMAX(k)/_dy) as Int);
      for (j in nrl..nrh) {
        for (i in ncl..nch) {
          mask(j,i) = 10000;
        }
      }
    }

    // 1.3 Re-check boundaries
    for (j in 1.._ncy) { // Left
      if (mask(j,1) != FLUID)
      mask(j,0) = 10000;
    }
    for (j in 1.._ncy) { // Right
      if (mask(j,_ncx) != FLUID)
      mask(j,_ncx+1) = 10000;
    }
    for (i in 0..(_ncx+1)) { // Top
      if (mask(_ncy,i) != FLUID)
      mask(_ncy+1,i) = 10000;
    }
    for (i in 0..(_ncx+1)) { // Bottom
      if (mask(1,i) != FLUID)
      mask(0,i) = 10000;
    }

    //2. set up inner domain non-fluid cell types
    for (j in 1.._ncy) {
      for (i in 1.._ncx) {
        if (mask(j,i) != FLUID) {
          if (mask(j,i+1) != FLUID) //East
          mask(j,i) += 1000;

          if (mask(j,i-1) != FLUID) //West
          mask(j,i) += 100;

          if (mask(j+1,i) != FLUID) //North
          mask(j,i) += 10;

          if (mask(j-1,i) != FLUID) //South
          mask(j,i) += 1;
        }
      }
    }
    return mask;
  }//end function

}//end class
