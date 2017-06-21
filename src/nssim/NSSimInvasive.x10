package nssim;
import x10.io.File;
import x10.io.FileWriter;
import x10.io.FileReader;
import x10.io.Printer;
import x10.compiler.Inline;
import x10.util.List;
import x10.util.ArrayList;

import invadeX10.invasic.Claim;
import invadeX10.invasic.IncarnationID;
import invadeX10.constraints.PEQuantity;
import invadeX10.constraints.ScalabilityHint;
import invadeX10.invasic.ProcessingElement;

public class NSSimInvasive extends NSSim{

  protected var claim:Claim;
  protected var world: PlaceGroup;
  protected var dist: Dist;
  protected var b: PlaceLocalHandle[NSBlock];

  /**********************
  *  Constructor
  * *********************/
  public def this()
  {
    super();
  }//end constructor

  public def run(): void
  {
    var dt: Double = 0.0;
    var t:  Double = 0.0;
    val t_end = final_time;
    var time_step: Int = 0;

    var vtk_freq:Double = 0.1;
    val vtk_cnt = GlobalCell.make[Int](1 as Int);
    val vtk_outfile = GlobalCell.make[String](output_name);

    // Create global mask and global fields
    val global_mask = GlobalRef[Array[Int]{rank==2}](create_geometry_mask());
    val global_u = GlobalRef[Array[Double]{rank==2}](new Array[Double]
    (0..((_ncy+1)) * 0..((_ncx+1)), _initial_velocity_x));
    val global_v = GlobalRef[Array[Double]{rank==2}](new Array[Double]
    (0..((_ncy+1)) * 0..((_ncx+1)), _initial_velocity_y));
    val global_p = GlobalRef[Array[Double]{rank==2}](new Array[Double]
    (0..((_ncy+1)) * 0..((_ncx+1)), _initial_pressure));
    val global_f = GlobalRef[Array[Double]{rank==2}](new Array[Double]
    (0..((_ncy+1)) * 0..((_ncx+1)), 0.0));
    val global_g = GlobalRef[Array[Double]{rank==2}](new Array[Double]
    (0..((_ncy+1)) * 0..((_ncx+1)), 0.0));

    ////////////////////////////////////////////////////////////////////////////
    // initial invasive operations
    ////////////////////////////////////////////////////////////////////////////
    //allocate the resources and initialize the claim
    var temp_map: Array[Int]{rank==1} = handle_resources_statically();
    val ID = temp_map;

    // create blocks of simulation to tun on the allocated resources (claim)
    b = PlaceLocalHandle.make[NSBlock](dist,
    ()=>create_block(ID, here.id, world.numPlaces(), _ncx, _ncy,
     global_mask, global_u, global_v, global_p, global_f, global_g));

    // Initialize U,V,P,F,G
    finish for (p in world) async at (p) {
      b().update_domain_uv();
      b().update_domain_p();
      b().update_domain_fg();
    }

    finish for (p in world) async at (p) {
      b().update_ghost_uv(b);
      b().update_ghost_p(b);
      b().update_ghost_fg(b);
    }

    // Initial VTK
    if (vtk_flag==1) {
      finish for (p in world) async at (p) {
        b().write_vts_file(vtk_outfile(), ID, 0, world.numPlaces());
      }
      b().write_vtm_file(vtk_outfile(), 0, world.numPlaces());
    }

    // Main loop
    while (t < t_end) {

      // 0 Invasive operations
      if (simulation_type==INVASIVE) {

        // 0.1 update everything in place 0 before all the invasion operations
        update_global_fields(b, claim, global_u, global_v, global_p, global_f,
        global_g);

        // 0.2 invasive operations
        temp_map = handle_resources_dynamically(time_step);
        val map =temp_map;

        // 0.3 create blocks of simulation to tun on the new allocated resources (claim)
        b = PlaceLocalHandle.make[NSBlock](dist,
        ()=>create_block(map, here.id, world.numPlaces(), _ncx, _ncy,
        global_mask, global_u, global_v, global_p, global_f, global_g));

      }

      // 1. compute dt
      finish for (p in world) async at (p) {
        b().compute_local_dt();
      }
      set_global_min_dt(b,claim);

      // 2. compoute F,G
      finish for (p in world) async at (p) {
        b().compute_fg();
        b().update_domain_fg();
      }
      finish for (p in world) async at (p) {
        b().update_ghost_fg(b);
      }

      // 3. solve for pressure P
      finish for (p in world) async at (p) {
        b().compute_p(b);
        b().update_domain_p();
      }
      finish for (p in world) async at (p) {
        b().update_ghost_p(b);
      }

      // 4. Compute velocity U,V
      finish for (p in world) async at (p) {
        b().compute_uv();
        b().update_domain_uv();
      }
      finish for (p in world) async at (p) {
        b().update_ghost_uv(b);
      }

      // 5. Update simulation time
      t += b().get_dt();
      Console.OUT.printf("Simulation completed at t = %.5f, dt = %.5f \n", t,
      b().get_dt());

      // 6. Write VTK output
      if (vtk_flag==1) {
        if ((Math.floor(t/vtk_freq) as Int) >= vtk_cnt()) {
          Console.OUT.printf("Writing output files at t = %.5f\n", t);

          val map = temp_map;
          finish for (p in world) async at (p) {
            b().write_vts_file(vtk_outfile(), map, vtk_cnt(), world.numPlaces());
          }
          b().write_vtm_file(vtk_outfile(), vtk_cnt(), world.numPlaces());
          vtk_cnt.set(vtk_cnt()+1);
        }
      }

      time_step+=1;
      Console.OUT.printf("at the end of time step %d, num proc = %d \n",
      time_step, world.numPlaces());
    }//end while

    claim.retreat();
    return;
  }// end funcion run()


  /********************
  *  Internal Methods
  ********************/

  protected def set_global_min_dt(b: PlaceLocalHandle[NSBlock], claim:Claim): void
  {
    val world: PlaceGroup = places(claim);
    val gdt = GlobalCell.make[Double](at(Place(0)) b().get_dt());
    var tmp: Double = 0.0;
    for (p in world) {
      tmp = at(p) b().get_dt();
      if (tmp < gdt()) gdt.set(tmp);
    }
    finish for (p in world) async at(p) {
      b().set_dt(gdt());
    }
    return;
  }//end function

  protected def update_global_fields(b: PlaceLocalHandle[NSBlock], claim:Claim,
  global_u: GlobalRef[Array[Double]{rank==2}],
  global_v: GlobalRef[Array[Double]{rank==2}],
  global_p: GlobalRef[Array[Double]{rank==2}],
  global_f: GlobalRef[Array[Double]{rank==2}],
  global_g: GlobalRef[Array[Double]{rank==2}]): void
  {
    var buffer: Array[Double]{rank==2};
    var ncx: Int;
    var ncy: Int;
    var gimin: Int;
    var gimax: Int;
    var gjmin: Int;
    var gjmax: Int;
    val world: PlaceGroup = places(claim);

    //update global u
    for (p in world) {

      gimin    = at (p) b().get_gimin();
      gimax    = at (p) b().get_gimax();
      gjmin    = at (p) b().get_gjmin();
      gjmax    = at (p) b().get_gjmax();
      ncx = gimax - gimin + 1;
      ncy = gjmax - gjmin + 1;

      // u
      buffer = at (p) b().get_u();
      for (gj in (gjmin-1)..(gjmax+1)) {
        for (gi in (gimin-1)..(gimax+1)) {
          global_u()(gj,gi) = buffer(lj(gj,gjmin),li(gi,gimin));
        }
      }
      // v
      buffer = at (p) b().get_v();
      for (gj in (gjmin-1)..(gjmax+1)) {
        for (gi in (gimin-1)..(gimax+1)) {
          global_v()(gj,gi) = buffer(lj(gj,gjmin),li(gi,gimin));
        }
      }
      // p
      buffer = at (p) b().get_p();
      for (gj in (gjmin-1)..(gjmax+1)) {
        for (gi in (gimin-1)..(gimax+1)) {
          global_p()(gj,gi) = buffer(lj(gj,gjmin),li(gi,gimin));
        }
      }
      // f
      buffer = at (p) b().get_f();
      for (gj in (gjmin-1)..(gjmax+1)) {
        for (gi in (gimin-1)..(gimax+1)) {
          global_f()(gj,gi) = buffer(lj(gj,gjmin),li(gi,gimin));
        }
      }
      // g
      buffer = at (p) b().get_g();
      for (gj in (gjmin-1)..(gjmax+1)) {
        for (gi in (gimin-1)..(gimax+1)) {
          global_g()(gj,gi) = buffer(lj(gj,gjmin),li(gi,gimin));
        }
      }

    }

    return;
  }//end function

  protected def places(claim:Claim):PlaceGroup{
    val pes = claim.processingElements();
    val cmp = (pe1:ProcessingElement, pe2:ProcessingElement) => {
      pe1.address.place().id.compareTo(pe2.address.place().id)
    };
    pes.sort(cmp);

    val list = new ArrayList[Place]();
    list.add(here);
    for(pe in pes) {
      val place = pe.address.place();
      if(!list.contains(place)) {
        list.add(place);
      }
    }
    return new SparsePlaceGroup(list.toArray().sequence());
  }

  // invesive resource managing methods
  // this function should be called to initialize the claim
  // and the other invasive classes
  def handle_resources_statically():Array[Int]{rank==1}{
    // NP is the bumber of processors set in Parameters
    val map = new Array[Int]( ( 0..(NP-1) ), 0);
    for ( i in (0..(NP-1)) ){
      map (i) = i;
    }
    claim = Claim.invade(new PEQuantity(NP-1));
    world = places(claim);
    dist = Dist.makeUnique(world);

    return map;
  }

  def handle_resources_dynamically(time_step:Int):Array[Int]{rank==1}{
    // processing number pool -> it is used as a
    // reference to manage the number of the processors
    // pool = {1, 2, 4, 8, 16}
    // non- invasive parallel simulation
    var pool: Array[Int] = new Array[Int](5 , (i:Int) => Math.pow2(i));
    var temp1: Int;
    // NP is the bumber of processors set in Parameters

    // 1 managing the processing count -> later on can be
    // done using a smart processing manager in each time
    // step
    temp1=NP;
    NP = pool(time_step%5);

    // 2 invasion process
    // invade in case more porcessors are needed
    if (NP>temp1){
      claim.reinvade(new PEQuantity(NP-temp1));
    }

    // retreat in case less porcessors are needed
    if (NP<temp1){
      for (i in 1..(temp1-NP)){
        claim.retreat(claim.processingElements()(claim.size()-1));
      }
    }

    // 3 invade
    world = places(claim);
    dist = Dist.makeUnique(world);

    // 4 "map" used to map the processor ids in the
    // in the simulation according to the ids
    // specified in the claim
    val map = new Array[Int]( ( 0..(world.numPlaces()-1) ), 0);
    var index:Int =0;
    for ( p in world ){
      map(index) = at (p) here.id;
      index++;
    }

    return map;
  }

  @Inline
  private def li(gi: Int, gimin: Int): Int {
    return gi - gimin + 1 as Int;
  }

  @Inline
  private def lj(gj: Int, gjmin: Int): Int {
    return gj - gjmin + 1 as Int;
  }

}//end class
