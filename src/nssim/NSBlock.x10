package nssim;
import x10.io.File;
import x10.io.FileWriter;
import x10.io.Printer;
import x10.compiler.Inline;

public class NSBlock extends Parameters {

  /**
  *=====================
  * Row-Major ordering:
  *=====================
  * 	  (y-dim)
  * 		M|           j: index over y-dim, outer loop
  * 	  j	 |           i: index over x-dim, inner loop
  * 		0|______
  *        0    N (x-dim)
  *           i
  *
  *    domain = new Array[Int]((0..M)*(0..N), 0);
  *
  * 	  domain(j,i) is the point at row-j, column-i
  *
  * NOTE:
  *  - Row-Major scheme is used for all rank==2 Arrays
  *  - Row-Major scheme is applied on domain decomposition (block ordering)
  *
  *
  * ======================================
  *  Index conversion based on Row-Major:
  * ======================================
  *
  *    seq = j*(N+1) + i   (See figure above)
  *    i = seq % (N+1)
  *    j = (seq / (N+1)) as Int
  *
  * ===============================
  *  Global | Local index mapping:
  * ===============================
  *
  *    Inner domain: [1, ncx  ] x [1, ncy  ], size = ncx * ncy
  *    Whole domain: [0, ncx+1] x [0, ncy+1], size = (ncx+2) * (ncy+2)
  *
  * E.g.: Simulation doamin 8x8: gncx=8, gncy=8
  *
  * 	  global domain: inner=[1,8]x[1,8], whole=[0,9]x[0,9]
  *
  *    local doamin assignment:
  *		gimin=5, gimax=7  =>  ncx = gimax-gimin+1 = 3
  * 		gjmin=3, gjmax=5  =>  ncy = gjmax-gjmin+1 = 3
  *
  * 	  local doamin:  inner=[1,3]x[1,3], whole=[0,4]x[0,4]
  *
  *      =>  gi = i + (gimin-1)   =>  i = gi - gimin + 1
  *          gj = j + (gjmin-1)       j = gj - gjmin + 1
  *
  *    Note: Global index prefixed with g: gi, gj, gncx, gncy.
  *          Local index no prefix: i, j, ncx, ncy
  */

  private val _id: Int;   // Block ID
  private val _nid: Int;	// North neighbor block ID
  private val _sid: Int;  // South neighbor block ID
  private val _eid: Int;  // East neighbor block ID
  private val _wid: Int;  // West neighbor block ID

  private val _gimin: Int; // Glboal index of local inner domain imin (x-dimension, inner loop)
  private val _gimax: Int; // Glboal index of local inner domain imax
  private val _gjmin: Int; // Glboal index of local inner domain jmin (y-dimension, outer loop)
  private val _gjmax: Int; // Glboal index of local inner domain jmax

  private val _M: Array[Int]{rank==2};	// Geometry masks
  private val _U: Array[Double]{rank==2}; // Velocity field in x-dimension
  private val _V: Array[Double]{rank==2}; // Velocity field in y-dimension
  private val _P: Array[Double]{rank==2}; // Pressure field
  private val _F: Array[Double]{rank==2};	// Intermediate step for U
  private val _G: Array[Double]{rank==2};	// Intermediate step for V

  private var _dt: Double;  // simulation time step

  /**********************
  *  Constructor overloaded
  * *********************/
  public def this(
  id: Int,
  nid: Int,
  sid: Int,
  eid: Int,
  wid: Int,
  gimin: Int,
  gimax: Int,
  gjmin: Int,
  gjmax: Int,
  gncx: Int,
  gncy: Int,
  global_mask: GlobalRef[Array[Int]{rank==2}],
  global_u: GlobalRef[Array[Double]{rank==2}],
  global_v: GlobalRef[Array[Double]{rank==2}],
  global_p: GlobalRef[Array[Double]{rank==2}],
  global_f: GlobalRef[Array[Double]{rank==2}],
  global_g: GlobalRef[Array[Double]{rank==2}])
  {

    super ();

    _ncx = gimax - gimin + 1;
    _ncy = gjmax - gjmin + 1;
    _dx = _domain_size_x / gncx;
    _dy = _domain_size_y / gncy;

    _id = id;
    _nid = (nid >= 0) ? nid : _boundary_north;
    _sid = (sid >= 0) ? sid : _boundary_south;
    _eid = (eid >= 0) ? eid : _boundary_east;
    _wid = (wid >= 0) ? wid : _boundary_west;
    _gimin = gimin;
    _gimax = gimax;
    _gjmin = gjmin;
    _gjmax = gjmax;

    _dt = 0.0;

    _M = new Array[Int](0..(_ncy+1) * 0..(_ncx+1), 0);
    for (gj in (gjmin-1)..(gjmax+1)) {
      for (gi in (gimin-1)..(gimax+1)) {
        _M(_lj(gj),_li(gi)) = at(global_mask) global_mask()(gj,gi);
      }
    }

    _U = new Array[Double](0..(_ncy+1) * 0..(_ncx+1), 0.0);
    for (gj in (gjmin-1)..(gjmax+1)) {
      for (gi in (gimin-1)..(gimax+1)) {
        _U(_lj(gj),_li(gi)) = at(global_u) global_u()(gj,gi);
      }
    }
    _V = new Array[Double](0..(_ncy+1) * 0..(_ncx+1), 0.0);
    for (gj in (gjmin-1)..(gjmax+1)) {
      for (gi in (gimin-1)..(gimax+1)) {
        _V(_lj(gj),_li(gi)) = at(global_v) global_v()(gj,gi);
      }
    }
    _P = new Array[Double](0..(_ncy+1) * 0..(_ncx+1), 0.0);
    for (gj in (gjmin-1)..(gjmax+1)) {
      for (gi in (gimin-1)..(gimax+1)) {
        _P(_lj(gj),_li(gi)) = at(global_p) global_p()(gj,gi);
      }
    }
    _F = new Array[Double](0..(_ncy+1) * 0..(_ncx+1), 0.0);
    for (gj in (gjmin-1)..(gjmax+1)) {
      for (gi in (gimin-1)..(gimax+1)) {
        _F(_lj(gj),_li(gi)) = at(global_f) global_f()(gj,gi);
      }
    }
    _G = new Array[Double](0..(_ncy+1) * 0..(_ncx+1), 0.0);
    for (gj in (gjmin-1)..(gjmax+1)) {
      for (gi in (gimin-1)..(gimax+1)) {
        _G(_lj(gj),_li(gi)) = at(global_g) global_g()(gj,gi);
      }
    }

    return;
  }//end constructor

  /**********************
  *  Computation dt
  **********************/
  public def set_dt(dt: Double): void
  {
    _dt = dt;
  }

  public def get_dt(): Double
  {
    return _dt;
  }

  public def compute_local_dt(): void
  {
    // compute local dt
    val abs_umax = get_matrix_fabs_max(_U, 1, _ncy, 0, _ncx);
    val abs_vmax = get_matrix_fabs_max(_V, 0, _ncy, 1, _ncx);
    _dt = _tau * Math.min(_re* _dx*_dx * _dy*_dy / (2.0 * (_dx*_dx + _dy*_dy)), Math.min(_dx/abs_umax, _dy/abs_vmax));
  }

  /**********************
  *  Computation F,G
  **********************/
  public def compute_fg(): void
  {
    //compute F
    for (j in 1.._ncy)
    for (i in 1.._ncx)
    _F(j,i) = _U(j,i) + _dt * ( laplacian(_U,j,i,_dx,_dy)/_re
    - FD_x_U2(_U,j,i,_dx,_alpha)
    - FD_y_UV(_U,_V,j,i,_dy,_alpha)
    + _external_force_x );

    // only on domain west/east boundaries
    if (_wid < 0) {
      for (j in 1.._ncy)
      _F(j,0) = _U(j,0);
    }
    if (_eid < 0) {
      for (j in 1.._ncy)
      _F(j,_ncx) = _U(j,_ncx);
    }

    //compute G
    for (j in 1.._ncy)
    for (i in 1.._ncx)
    _G(j,i) = _V(j,i) + _dt * ( laplacian(_V,j,i,_dx,_dy)/_re
    - FD_x_UV(_U,_V,j,i,_dx,_alpha)
    - FD_y_V2(_V,j,i,_dy,_alpha)
    + _external_force_y );

    // only on domain south/north boundaries
    if (_sid < 0) {
      for (i in 1.._ncx)
      _G(0,i) = _V(0,i);
    }
    if (_nid < 0) {
      for (i in 1.._ncx)
      _G(_ncy,i) = _V(_ncy,i);
    }
    return;
  }//end function

  public def update_domain_fg(): void
  {
    for (j in 1.._ncy) {
      for (i in 1.._ncx) {
        if (_M(j,i) == FLUID) {
          continue;
        } else if (_M(j,i) == B_N) {
          _G(j,i) = _V(j,i);
        } else if (_M(j,i) == B_S) {
          _G(j,i) = _V(j,i);
          _G(j-1,i) = _V(j-1,i);
        } else if (_M(j,i) == B_E) {
          _F(j,i) = _U(j,i);
        } else if (_M(j,i) == B_W) {
          _F(j,i) = _U(j,i);
          _F(j,i-1) = _U(j,i-1);
        } else if (_M(j,i) == B_NE) {
          _F(j,i) = _U(j,i);
          _G(j,i) = _V(j,i);
        } else if (_M(j,i) == B_SE) {
          _F(j,i) = _U(j,i);
          _G(j-1,i) = _V(j-1,i);
        } else if (_M(j,i) == B_NW) {
          _F(j,i-1) = _U(j,i-1);
          _G(j,i) = _V(j,i);
        } else if (_M(j,i) == B_SW) {
          _F(j,i-1) = _U(j,i-1);
          _G(j-1,i) = _V(j-1,i);
        } else if (_M(j,i) == B_IN) {
          _F(j,i) = _U(j,i);
          _G(j,i) = _V(j,i);
        } else { //forbidden cases
          Console.OUT.printf("Forbidden cell. Probably caused by low resolution. Simulation abort!"+"\n update_obstacles_fg: i=%d, j=%d, _M=%d \n",i,j,_M(j,i));
        }
      }
    }
    return;
  }//end function

  public def update_ghost_fg(neiblock: PlaceLocalHandle[NSBlock]): void
  {
    //north
    if (_nid >= 0) {
      val fbuf = at(Place(_nid)) neiblock().get_f_south();
      for (i in 0.._ncx)
      _F(_ncy+1,i) = fbuf(i);

      val gbuf = at(Place(_nid)) neiblock().get_g_south();
      for (i in 1.._ncx)
      _G(_ncy+1,i) = gbuf(i);
    }
    //south
    if (_sid >= 0) {
      val fbuf = at(Place(_sid)) neiblock().get_f_north();
      for (i in 0.._ncx)
      _F(0,i) = fbuf(i);

      val gbuf = at(Place(_sid)) neiblock().get_g_north();
      for (i in 1.._ncx)
      _G(0,i) = gbuf(i);
    }
    //east
    if (_eid >= 0) {
      val fbuf = at(Place(_eid)) neiblock().get_f_west();
      for (j in 1.._ncy)
      _F(j,_ncx+1) = fbuf(j);

      val gbuf = at(Place(_eid)) neiblock().get_g_west();
      for (j in 0.._ncy)
      _G(j,_ncx+1) = gbuf(j);
    }
    //west
    if (_wid >= 0) {
      val fbuf = at(Place(_wid)) neiblock().get_f_east();
      for (j in 1.._ncy)
      _F(j,0) = fbuf(j);

      val gbuf = at(Place(_wid)) neiblock().get_g_east();
      for (j in 0.._ncy)
      _G(j,0) = gbuf(j);
    }
  }//end function


  /**********************
  *  Computation P
  * *********************/
  //parallel solve
  public def compute_p(neiblock: PlaceLocalHandle[NSBlock]): void
  {
    solve(_ncx, _ncy, _dx, _dy, _dt, _omega, _F, _G, _P, neiblock);
  }//end function

  //serial solve
  public def compute_p(): void
  {
    solve(_ncx, _ncy, _dx, _dy, _dt, _omega, _F, _G, _P);
  }//end function

  public def update_domain_p(): void
  {
    for (j in 1.._ncy) {
      for (i in 1.._ncx) {
        if(_M(j,i) == FLUID) {
          continue;
        } else if (_M(j,i) == B_N) {
          _P(j,i) = _P(j+1,i);
        } else if (_M(j,i) == B_S) {
          _P(j,i) = _P(j-1,i);
        } else if (_M(j,i) == B_E) {
          _P(j,i) = _P(j,i+1);
        } else if (_M(j,i) == B_W) {
          _P(j,i) = _P(j,i-1);
        } else if (_M(j,i) == B_NE) {
          _P(j,i) = (_P(j+1,i) + _P(j,i+1)) / 2.0;
        } else if (_M(j,i) == B_SE) {
          _P(j,i) = (_P(j-1,i) + _P(j,i+1)) / 2.0;
        } else if (_M(j,i) == B_NW) {
          _P(j,i) = (_P(j+1,i) + _P(j,i-1)) / 2.0;
        } else if (_M(j,i) == B_SW) {
          _P(j,i) = (_P(j-1,i) + _P(j,i-1)) / 2.0;
        } else if (_M(j,i) == B_IN) {
          _P(j,i) = 0;
        } else { //forbidden cases
          Console.OUT.printf("Forbidden cell. Probably caused by low resolution. Simulation abort! "+"\n update_obstacles_p: i=%d, j=%d, _M=%d \n",i,j,_M(j,i));
        }
      }
    }
    return;
  }//end function

  public def update_ghost_p(neiblock: PlaceLocalHandle[NSBlock]): void
  {
    //north
    if (_nid >= 0) {
      val pbuf = at(Place(_nid)) neiblock().get_p_south();
      for (i in 1.._ncx)
      _P(_ncy+1,i) = pbuf(i);
    }
    //south
    if (_sid >= 0) {
      val pbuf = at(Place(_sid)) neiblock().get_p_north();
      for (i in 1.._ncx)
      _P(0,i) = pbuf(i);
    }
    //east
    if (_eid >= 0) {
      val pbuf = at(Place(_eid)) neiblock().get_p_west();
      for (j in 1.._ncy)
      _P(j,_ncx+1) = pbuf(j);
    }
    //west
    if (_wid >= 0) {
      val pbuf = at(Place(_wid)) neiblock().get_p_east();
      for (j in 1.._ncy)
      _P(j,0) = pbuf(j);
    }
  }//end function

  // SOR solver for pressure
  // parallel
  public def solve(
  ncx: Int,
  ncy: Int,
  dx: Double,
  dy: Double,
  dt: Double,
  omega: Double,
  F: Array[Double]{rank==2},
  G: Array[Double]{rank==2},
  P: Array[Double]{rank==2},
  neiblock: PlaceLocalHandle[NSBlock]): void
  {
    val ITERMAX = 10000 as Int;
    val TOLERANCE = 0.0001 as Double;
    val inv_dx2 = (1.0 / (dx * dx)) as Double;
    val inv_dy2 = (1.0 / (dy * dy)) as Double;
    val inv_dt_dx = (1.0 / dt / dx) as Double;
    val inv_dt_dy = (1.0 / dt / dy) as Double;
    val a = (1.0 - omega) as Double;
    val b = (0.5 * omega / (inv_dx2 + inv_dy2)) as Double;

    var res: Double;
    var tmp: Double;
    for (it in 1..ITERMAX) {
      // Swip over P
      for (j in 1..ncy) {
        for (i in 1..ncx) {
          P(j,i) = a * P(j,i) + b * ( inv_dx2*(P(j,i+1)+P(j,i-1))
          + inv_dy2*(P(j+1,i)+P(j-1,i))
          - (inv_dt_dx*(F(j,i)-F(j,i-1)) + inv_dt_dy*(G(j,i)-G(j-1,i))) );
        }
      }
      update_ghost_p(neiblock);

      // Compute residual
      res = 0.0;
      for (j in 1..ncy) {
        for (i in 1..ncx) {
          tmp = inv_dx2*(P(j,i+1)-2.0*P(j,i)+P(j,i-1)) +
          inv_dy2*(P(j+1,i)-2.0*P(j,i)+P(j-1,i)) -
          (inv_dt_dx*(F(j,i)-F(j,i-1)) + inv_dt_dy*(G(j,i)-G(j-1,i)));
          res += tmp*tmp;
        }
      }
      res = Math.sqrt(res/((ncy*ncx) as Double));
      if (res <= TOLERANCE) return;
    }
    Console.OUT.printf("SOR solver did not converge for %d iterations!\n", ITERMAX);
    return;
  }//end function

  // SOR solver for pressure
  // serial
  public def solve(
  ncx: Int,
  ncy: Int,
  dx: Double,
  dy: Double,
  dt: Double,
  omega: Double,
  F: Array[Double]{rank==2},
  G: Array[Double]{rank==2},
  P: Array[Double]{rank==2}): void
  {
    val ITERMAX = 10000 as Int;
    val TOLERANCE = 0.0001 as Double;
    val inv_dx2 = (1.0 / (dx * dx)) as Double;
    val inv_dy2 = (1.0 / (dy * dy)) as Double;
    val inv_dt_dx = (1.0 / dt / dx) as Double;
    val inv_dt_dy = (1.0 / dt / dy) as Double;
    val a = (1.0 - omega) as Double;
    val b = (0.5 * omega / (inv_dx2 + inv_dy2)) as Double;

    var res: Double;
    var tmp: Double;
    for (it in 1..ITERMAX) {
      // Swip over P
      for (j in 1..ncy) {
        for (i in 1..ncx) {
          P(j,i) = a * P(j,i) + b * ( inv_dx2*(P(j,i+1)+P(j,i-1))
          + inv_dy2*(P(j+1,i)+P(j-1,i))
          - (inv_dt_dx*(F(j,i)-F(j,i-1)) + inv_dt_dy*(G(j,i)-G(j-1,i))) );
        }
      }

      // Compute residual
      res = 0.0;
      for (j in 1..ncy) {
        for (i in 1..ncx) {
          tmp = inv_dx2*(P(j,i+1)-2.0*P(j,i)+P(j,i-1)) +
          inv_dy2*(P(j+1,i)-2.0*P(j,i)+P(j-1,i)) -
          (inv_dt_dx*(F(j,i)-F(j,i-1)) + inv_dt_dy*(G(j,i)-G(j-1,i)));
          res += tmp*tmp;
        }
      }
      res = Math.sqrt(res/((ncy*ncx) as Double));
      if (res <= TOLERANCE) return;
    }
    Console.OUT.printf("SOR solver did not converge for %d iterations!\n", ITERMAX);
    return;
  }//end function

  /**********************
  *  Computation U,V
  * *********************/
  public def compute_uv(): void
  {
    for (j in 1.._ncy)
    for (i in 1.._ncx)
    _U(j,i) = _F(j,i) - _dt * (_P(j,i+1) - _P(j,i)) / _dx;

    for (j in 1.._ncy)
    for (i in 1.._ncx)
    _V(j,i) = _G(j,i) - _dt * (_P(j+1,i) - _P(j,i)) / _dy;
    return;
  }//end function

  public def update_domain_uv(): void
  {
    for (j in 1.._ncy) {
      for (i in 1.._ncx) {
        if(_M(j,i) == FLUID) {
          continue;
        } else if (_M(j,i) == B_N) { // North edge cell
          _U(j,i) = -_U(j+1,i);
          _V(j,i) = 0.0;
        } else if (_M(j,i) == B_S) { // South edge cell
          _U(j,i) = -_U(j-1,i);
          _V(j,i) = 0.0;
          _V(j-1,i) = 0.0;
        } else if (_M(j,i) == B_E) { // East edge cell
          _U(j,i) = 0.0;
          _V(j,i) = -_V(j,i+1);
        } else if (_M(j,i) == B_W) { // West edge cell
          _U(j,i) = 0.0;
          _V(j,i) = -_V(j,i-1);
          _U(j,i-1) = 0.0;
        } else if (_M(j,i) == B_NE) { // North-east corner cell
          _U(j,i) = 0.0;
          _V(j,i) = 0.0;
        } else if (_M(j,i) == B_SE) { // South-east corner cell
          _U(j,i) = 0.0;
          _V(j,i) = -_V(j,i+1);
          _V(j-1,i) = 0.0;
        } else if (_M(j,i) == B_NW) { // North-west corner cell
          _U(j,i) = -_U(j+1,i);
          _V(j,i) = 0.0;
          _U(j,i-1) = 0.0;
        } else if (_M(j,i) == B_SW) { // South-west corner cell
          _U(j,i) = -_U(j-1,i);
          _V(j,i) = -_V(j,i-1);
          _U(j,i-1) = 0.0;
          _V(j-1,i) = 0.0;
        } else if (_M(j,i) == B_IN) {
          _U(j,i) = 0.0;
          _V(j,i) = 0.0;
        } else { //forbidden cases
          Console.OUT.printf("Forbidden cell. Probably caused by low resolution. Simulation abort!"+"\n update_obstacles_uv: i=%d, j=%d, _M=%d \n",i,j,_M(j,i));
        }
      }
    }
    return;
  }//end function

  //parallel
  public def update_ghost_uv(neiblock: PlaceLocalHandle[NSBlock]): void
  {
    //north
    if (_nid < 0) { //if north is domain boundary
      if (_boundary_north == BOUNDARY_TYPE_INLET) {
        for (i in 0.._ncx)
        _U(_ncy+1,i) = _inlet_velocity_x*2.0 - _U(_ncy,i);
        for (i in 1.._ncx)
        _V(_ncy,i) = _inlet_velocity_y;
      }
      else if (_boundary_north == BOUNDARY_TYPE_OUTLET) {
        for (i in 0.._ncx)
        _U(_ncy+1,i) = _U(_ncy,i);
        for (i in 1.._ncx)
        _V(_ncy,i) = _V(_ncy-1,i);
      }
      else if (_boundary_north == BOUNDARY_TYPE_NOSLIP) {
        for (i in 0.._ncx)
        _U(_ncy+1,i) = -_U(_ncy,i);
        for (i in 1.._ncx)
        _V(_ncy,i) = 0.0;
      }
      else if (_boundary_north == BOUNDARY_TYPE_FREESLIP) {
        for (i in 0.._ncx)
        _U(_ncy+1,i) = _U(_ncy,i);
        for (i in 1.._ncx)
        _V(_ncy,i) = 0.0;
      }
    } else {
      val ubuf = at(Place(_nid)) neiblock().get_u_south();
      for (i in 0.._ncx)
      _U(_ncy+1,i) = ubuf(i);

      val vbuf = at(Place(_nid)) neiblock().get_v_south();
      for (i in 1.._ncx)
      _V(_ncy+1,i) = vbuf(i);
    }

    //south
    if (_sid < 0) { //if south is domain boundary
      if (_boundary_south == BOUNDARY_TYPE_INLET) {
        for (i in 0.._ncx)
        _U(0,i) = _inlet_velocity_x*2.0 - _U(1,i);
        for (i in 1.._ncx)
        _V(0,i) = _inlet_velocity_y;
      }
      else if (_boundary_south == BOUNDARY_TYPE_OUTLET) {
        for (i in 0.._ncx)
        _U(0,i) = _U(1,i);
        for (i in 1.._ncx)
        _V(0,i) = _V(1,i);
      }
      else if (_boundary_south == BOUNDARY_TYPE_NOSLIP) {
        for (i in 0.._ncx)
        _U(0,i) = -_U(1,i);
        for (i in 1.._ncx)
        _V(0,i) = 0.0;
      }
      else if (_boundary_south == BOUNDARY_TYPE_FREESLIP) {
        for (i in 0.._ncx)
        _U(0,i) = _U(1,i);
        for (i in 1.._ncx)
        _V(0,i) = 0.0;
      }
    } else {
      val ubuf = at(Place(_sid)) neiblock().get_u_north();
      for (i in 0.._ncx)
      _U(0,i) = ubuf(i);

      val vbuf = at(Place(_sid)) neiblock().get_v_north();
      for (i in 1.._ncx)
      _V(0,i) = vbuf(i);
    }

    //east
    if (_eid < 0) { //if east is domain boundary
      if (_boundary_east == BOUNDARY_TYPE_INLET) {
        for (j in 1.._ncy)
        _U(j,_ncx) = _inlet_velocity_x;
        for (j in 0.._ncy)
        _V(j,_ncx+1) = _inlet_velocity_y*2.0 - _V(j,_ncx);
      }
      else if (_boundary_east == BOUNDARY_TYPE_OUTLET) {
        for (j in 1.._ncy)
        _U(j,_ncx) = _U(j,_ncx-1);
        for (j in 0.._ncy)
        _V(j,_ncx+1) = _V(j,_ncx);
      }
      else if (_boundary_east == BOUNDARY_TYPE_NOSLIP) {
        for (j in 1.._ncy)
        _U(j,_ncx) = 0.0;
        for (j in 0.._ncy)
        _V(j,_ncx+1) = -_V(j,_ncx);
      }
      else if (_boundary_east == BOUNDARY_TYPE_FREESLIP) {
        for (j in 1.._ncy)
        _U(j,_ncx) = 0.0;
        for (j in 0.._ncy)
        _V(j,_ncx+1) = _V(j,_ncx);
      }
    } else {
      val ubuf = at(Place(_eid)) neiblock().get_u_west();
      for (j in 1.._ncy)
      _U(j,_ncx+1) = ubuf(j);

      val vbuf = at(Place(_eid)) neiblock().get_v_west();
      for (j in 0.._ncy)
      _V(j,_ncx+1) = vbuf(j);
    }

    //west
    if (_wid < 0) { //if west is domain boundary
      if (_boundary_west == BOUNDARY_TYPE_INLET) {
        for (j in 1.._ncy)
        _U(j,0) = _inlet_velocity_x;
        for (j in 0.._ncy)
        _V(j,0) = _inlet_velocity_y*2.0 - _V(j,1);
      }
      else if (_boundary_west == BOUNDARY_TYPE_OUTLET) {
        for (j in 1.._ncy)
        _U(j,0) = _U(j,1);
        for (j in 0.._ncy)
        _V(j,0) = _V(j,1);
      }
      else if (_boundary_west == BOUNDARY_TYPE_NOSLIP) {
        for (j in 1.._ncy)
        _U(j,0) = 0.0;
        for (j in 0.._ncy)
        _V(j,0) = -_V(j,1);
      }
      else if (_boundary_west == BOUNDARY_TYPE_FREESLIP) {
        for (j in 1.._ncy)
        _U(j,0) = 0.0;
        for (j in 0.._ncy)
        _V(j,0) = _V(j,1);
      }
    } else {
      val ubuf = at(Place(_wid)) neiblock().get_u_east();
      for (j in 1.._ncy)
      _U(j,0) = ubuf(j);

      val vbuf = at(Place(_wid)) neiblock().get_v_east();
      for (j in 0.._ncy)
      _V(j,0) = vbuf(j);
    }
  }//end function

  //just for serial
  public def update_boundary_uv(): void
  {
    //north
    if (_boundary_north == BOUNDARY_TYPE_INLET) {
      for (i in 0.._ncx)
      _U(_ncy+1,i) = _inlet_velocity_x*2.0 - _U(_ncy,i);
      for (i in 1.._ncx)
      _V(_ncy,i) = _inlet_velocity_y;
    }
    else if (_boundary_north == BOUNDARY_TYPE_OUTLET) {
      for (i in 0.._ncx)
      _U(_ncy+1,i) = _U(_ncy,i);
      for (i in 1.._ncx)
      _V(_ncy,i) = _V(_ncy-1,i);
    }
    else if (_boundary_north == BOUNDARY_TYPE_NOSLIP) {
      for (i in 0.._ncx)
      _U(_ncy+1,i) = -_U(_ncy,i);
      for (i in 1.._ncx)
      _V(_ncy,i) = 0.0;
    }
    else if (_boundary_north == BOUNDARY_TYPE_FREESLIP) {
      for (i in 0.._ncx)
      _U(_ncy+1,i) = _U(_ncy,i);
      for (i in 1.._ncx)
      _V(_ncy,i) = 0.0;
    }

    //south
    if (_boundary_south == BOUNDARY_TYPE_INLET) {
      for (i in 0.._ncx)
      _U(0,i) = _inlet_velocity_x*2.0 - _U(1,i);
      for (i in 1.._ncx)
      _V(0,i) = _inlet_velocity_y;
    }
    else if (_boundary_south == BOUNDARY_TYPE_OUTLET) {
      for (i in 0.._ncx)
      _U(0,i) = _U(1,i);
      for (i in 1.._ncx)
      _V(0,i) = _V(1,i);
    }
    else if (_boundary_south == BOUNDARY_TYPE_NOSLIP) {
      for (i in 0.._ncx)
      _U(0,i) = -_U(1,i);
      for (i in 1.._ncx)
      _V(0,i) = 0.0;
    }
    else if (_boundary_south == BOUNDARY_TYPE_FREESLIP) {
      for (i in 0.._ncx)
      _U(0,i) = _U(1,i);
      for (i in 1.._ncx)
      _V(0,i) = 0.0;
    }

    //east
    if (_boundary_east == BOUNDARY_TYPE_INLET) {
      for (j in 1.._ncy)
      _U(j,_ncx) = _inlet_velocity_x;
      for (j in 0.._ncy)
      _V(j,_ncx+1) = _inlet_velocity_y*2.0 - _V(j,_ncx);
    }
    else if (_boundary_east == BOUNDARY_TYPE_OUTLET) {
      for (j in 1.._ncy)
      _U(j,_ncx) = _U(j,_ncx-1);
      for (j in 0.._ncy)
      _V(j,_ncx+1) = _V(j,_ncx);
    }
    else if (_boundary_east == BOUNDARY_TYPE_NOSLIP) {
      for (j in 1.._ncy)
      _U(j,_ncx) = 0.0;
      for (j in 0.._ncy)
      _V(j,_ncx+1) = -_V(j,_ncx);
    }
    else if (_boundary_east == BOUNDARY_TYPE_FREESLIP) {
      for (j in 1.._ncy)
      _U(j,_ncx) = 0.0;
      for (j in 0.._ncy)
      _V(j,_ncx+1) = _V(j,_ncx);
    }

    //west
    if (_boundary_west == BOUNDARY_TYPE_INLET) {
      for (j in 1.._ncy)
      _U(j,0) = _inlet_velocity_x;
      for (j in 0.._ncy)
      _V(j,0) = _inlet_velocity_y*2.0 - _V(j,1);
    }
    else if (_boundary_west == BOUNDARY_TYPE_OUTLET) {
      for (j in 1.._ncy)
      _U(j,0) = _U(j,1);
      for (j in 0.._ncy)
      _V(j,0) = _V(j,1);
    }
    else if (_boundary_west == BOUNDARY_TYPE_NOSLIP) {
      for (j in 1.._ncy)
      _U(j,0) = 0.0;
      for (j in 0.._ncy)
      _V(j,0) = -_V(j,1);
    }
    else if (_boundary_west == BOUNDARY_TYPE_FREESLIP) {
      for (j in 1.._ncy)
      _U(j,0) = 0.0;
      for (j in 0.._ncy)
      _V(j,0) = _V(j,1);
    }

  }//end function

  /**********************************
  *  Ghost Layer communication
  * ********************************/
  // Array _U:
  public def get_u_north(): Array[Double]{rank==1} {
    return new Array[Double](0.._ncx, (p:Point)=>{ _U(_ncy,p(0)) });
  }
  public def get_u_south(): Array[Double]{rank==1} {
    return new Array[Double](0.._ncx, (p:Point)=>{ _U(1,p(0)) });
  }
  public def get_u_east(): Array[Double]{rank==1} {
    return new Array[Double](1.._ncy, (p:Point)=>{ _U(p(0),_ncx) });
  }
  public def get_u_west(): Array[Double]{rank==1} {
    return new Array[Double](1.._ncy, (p:Point)=>{ _U(p(0),1) });
  }
  // Array _V:
  public def get_v_north(): Array[Double]{rank==1} {
    return new Array[Double](1.._ncx, (p:Point)=>{ _V(_ncy,p(0)) });
  }
  public def get_v_south(): Array[Double]{rank==1} {
    return new Array[Double](1.._ncx, (p:Point)=>{ _V(1,p(0)) });
  }
  public def get_v_east(): Array[Double]{rank==1} {
    return new Array[Double](0.._ncy, (p:Point)=>{ _V(p(0),_ncx) });
  }
  public def get_v_west(): Array[Double]{rank==1} {
    return new Array[Double](0.._ncy, (p:Point)=>{ _V(p(0),1) });
  }
  // Array _P:
  public def get_p_north(): Array[Double]{rank==1} {
    return new Array[Double](1.._ncx, (p:Point)=>{ _P(_ncy,p(0)) });
  }
  public def get_p_south(): Array[Double]{rank==1} {
    return new Array[Double](1.._ncx, (p:Point)=>{ _P(1,p(0)) });
  }
  public def get_p_east(): Array[Double]{rank==1} {
    return new Array[Double](1.._ncy, (p:Point)=>{ _P(p(0),_ncx) });
  }
  public def get_p_west(): Array[Double]{rank==1} {
    return new Array[Double](1.._ncy, (p:Point)=>{ _P(p(0),1) });
  }
  // Array _F:
  public def get_f_north(): Array[Double]{rank==1} {
    return new Array[Double](0.._ncx, (p:Point)=>{ _F(_ncy,p(0)) });
  }
  public def get_f_south(): Array[Double]{rank==1} {
    return new Array[Double](0.._ncx, (p:Point)=>{ _F(1,p(0)) });
  }
  public def get_f_east(): Array[Double]{rank==1} {
    return new Array[Double](1.._ncy, (p:Point)=>{ _F(p(0),_ncx) });
  }
  public def get_f_west(): Array[Double]{rank==1} {
    return new Array[Double](1.._ncy, (p:Point)=>{ _F(p(0),1) });
  }
  // Array _G:
  public def get_g_north(): Array[Double]{rank==1} {
    return new Array[Double](1.._ncx, (p:Point)=>{ _G(_ncy,p(0)) });
  }
  public def get_g_south(): Array[Double]{rank==1} {
    return new Array[Double](1.._ncx, (p:Point)=>{ _G(1,p(0)) });
  }
  public def get_g_east(): Array[Double]{rank==1} {
    return new Array[Double](0.._ncy, (p:Point)=>{ _G(p(0),_ncx) });
  }
  public def get_g_west(): Array[Double]{rank==1} {
    return new Array[Double](0.._ncy, (p:Point)=>{ _G(p(0),1) });
  }


  /**********************************
  *  Domain communication / update the global fields
  * ********************************/

  // U
  public def get_u(): Array[double]{rank==2} {
    //return new Array[double](0.._ncy * 0.._ncx, (p:Point)=>{ _U(p(0),p(1)) });
    return _U;
  }
  // V
  public def get_v(): Array[double]{rank==2} {
    //return new Array[double](0.._ncy * 0.._ncx, (p:Point)=>{ _V(p(0),p(1)) });
    return _V;
  }
  // P
  public def get_p(): Array[double]{rank==2} {
    //return new Array[double](0.._ncy * 0.._ncx, (p:Point)=>{ _P(p(0),p(1)) });
    return _P;
  }
  // F
  public def get_f(): Array[double]{rank==2} {
    //return new Array[double](0.._ncy * 0.._ncx, (p:Point)=>{ _F(p(0),p(1)) });
    return _F;
  }
  // G
  public def get_g(): Array[double]{rank==2} {
    //return new Array[double](0.._ncy * 0.._ncx, (p:Point)=>{ _G(p(0),p(1)) });
    return _G;
  }

  /**********************************
  *  Getters
  * ********************************/

  public def get_gimin(): Int {
    return this._gimin;
  }

  public def get_gimax(): Int {
    return this._gimax;
  }

  public def get_gjmin(): Int {
    return this._gjmin;
  }

  public def get_gjmax(): Int {
    return this._gjmax;
  }


  public def collect_at_zero(
  global_u: GlobalRef[Array[Double]{rank==2}],
  global_v: GlobalRef[Array[Double]{rank==2}],
  global_p: GlobalRef[Array[Double]{rank==2}],
  global_f: GlobalRef[Array[Double]{rank==2}],
  global_g: GlobalRef[Array[Double]{rank==2}])
  {
    for (gj in (_gjmin-1)..(_gjmax+1)) {
      for (gi in (_gimin-1)..(_gimax+1)) {
        at(global_u) global_u()(gj,gi) = _U(_lj(gj),_li(gi));
      }
    }
  }

  /**********************************
  *  Output writers, VTS and VTM
  * ********************************/
  public def write_vts_header_coord(printer: Printer): void
  {
    printer.printf("<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
    printer.printf("<StructuredGrid WholeExtent=\"0 %d 0 %d 0 0\">\n",_ncx,_ncy);
    printer.printf("<Piece Extent=\"0 %d 0 %d 0 0\">\n", _ncx, _ncy);

    printer.printf("<Points>\n");

    val originx = ((_gimin-1) * _dx) as Double;
    val originy = ((_gjmin-1) * _dy) as Double;

    printer.printf("<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n", originx, originy);
    for (j in 0.._ncy)
    for (i in 0.._ncx)
    printer.printf("%f %f %f\n", (originx + i*_dx), (originy + j*_dy), 0.0);
    printer.printf("</DataArray>\n");
    printer.printf("</Points>\n");
    return;
  }//end function

  // output writer
  public def write_vts_file(
  szProblem: String,
  map: Array[Int]{rank==1},
  timeStepNumber: Int,
  nprocs: Int): void
  {
    val outfile:String;

    //map inverse
    var rank: Int;  // rank
    rank=0;
    for (i in (0..(nprocs-1))  ){
      if (map(i)==here.id)
      rank = i;
    }

    outfile = szProblem+"_"+(timeStepNumber.toString()).trim()+"_"+rank+".vts";
    val file = new File(outfile);
    val writer = new FileWriter(file);
    val printer = new Printer(writer);

    //write header & coordinates
    write_vts_header_coord(printer);

    //write velocity
    printer.print("<PointData Vectors=\"velocity\">\n");
    printer.printf("<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">");
    for (j in 0.._ncy) {
      for (i in 0.._ncx) {
        if (is_point_in_obstacle_cell(_M,j,i))
        printer.printf("%f %f %f\n", 0.0, 0.0, 0.0);
        else
        printer.printf("%f %f %f\n", ((_U(j,i)+_U(j,i+1))*0.5), ((_V(j,i)+_V(j+1,i))*0.5), 0.0);
      }
    }
    printer.printf("</DataArray>\n");
    printer.printf("</PointData>\n");

    //write pressure
    printer.printf("<CellData Scalars=\"pressure\">\n");
    printer.printf("<DataArray type=\"Float32\" Name=\"pressure\" format=\"ascii\">\n");
    for (j in 1.._ncy) {
      for (i in 1.._ncx) {
        if (_M(j,i) == FLUID)
        printer.printf("%f\n", _P(j,i));
        else
        printer.printf("%f\n", 0.0);
      }
    }
    printer.printf("</DataArray>\n");
    printer.printf("</CellData>\n");
    printer.printf("</Piece>\n");
    printer.printf("</StructuredGrid>\n");
    printer.printf("</VTKFile>");

    printer.flush();
    writer.flush();
  }//end function

  public def write_vtm_file(
  szProblem: String,
  timeStepNumber: Int,
  nprocs: Int): void
  {
    val outfile:String;

    outfile = szProblem+"_"+(timeStepNumber.toString()).trim()+".vtm";
    val file = new File(outfile);
    val writer = new FileWriter(file);
    val printer = new Printer(writer);

    printer.printf("<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
    printer.printf("<vtkMultiBlockDataSet>\n");

    for (i in (0..(nprocs-1))  ){
      printer.printf("<DataSet index=\"%d\" file=\"artery_%s_%d.vts\">\n",i,(timeStepNumber.toString()).trim(),i);
      printer.printf("</DataSet>\n");
    }

    printer.printf("</vtkMultiBlockDataSet>\n");
    printer.printf("</VTKFile>");

    printer.flush();
    writer.flush();

  }

  protected def get_matrix_fabs_max(
  m: Array[Double]{rank==2},
  nrl: Int, nrh: Int, ncl: Int, nch: Int): Double
  {
    var newmax:Double;
    newmax = Math.abs(m(nrl,ncl));
    for (j in nrl..nrh) {
      for (i in ncl..nch) {
        if (newmax < Math.abs(m(j,i)))
        newmax = Math.abs(m(j,i));
      }
    }
    return newmax;
  }//end function

  /*****************
  *  DEBUG only
  *****************/
  public def print_block_info(): void
  {
    Console.OUT.printf("\nAt place %d: Block %d\n", here.id, _id);
    Console.OUT.printf("Neighbors: north %d, south %d, east %d, west %d\n", _nid, _sid, _eid, _wid);
    Console.OUT.printf("Index range: i in [%d, %d], j in [%d, %d]\n\n", _gimin, _gimax, _gjmin, _gjmax);
    //print_mask(_M);
  }//end function


  /**********************
  *  Inline Methods
  * *********************/

  @Inline
  protected def laplacian(
  m: Array[Double]{rank==2},
  j: Int,
  i: Int,
  dx: Double,
  dy: Double): Double
  {
    return (m(j,i+1) - 2.0*m(j,i) + m(j,i-1)) / (dx*dx) +
    (m(j+1,i) - 2.0*m(j,i) + m(j-1,i)) / (dy*dy);
  }

  @Inline
  protected def FD_x_U2(
  u: Array[Double]{rank==2},
  j: Int,
  i: Int,
  dx: Double,
  alpha: Double): Double
  {
    return (            ( (u(j,  i) + u(j,i+1)) * (u(j,  i) + u(j,i+1)) -
    (u(j,i-1) + u(j,  i)) * (u(j,i-1) + u(j,  i))
    )
    + alpha * ( Math.abs(u(j,  i) + u(j,i+1)) * (u(j,  i) - u(j,i+1)) -
    Math.abs(u(j,i-1) + u(j,  i)) * (u(j,i-1) - u(j,  i))
    )
    ) / (dx*4.0);
  }

  @Inline
  protected def FD_y_V2(
  v: Array[Double]{rank==2},
  j: Int,
  i: Int,
  dy: Double,
  alpha: Double): Double
  {
    return (           ( (v(j,i) + v(j+1,i)) * (v(j,i) + v(j+1,i)) -
    (v(j-1,i) + v(j,i)) * (v(j-1,i) + v(j,i))
    )
    + alpha * ( Math.abs(v(j  ,i) + v(j+1,i)) * (v(j  ,i) - v(j+1,i)) -
    Math.abs(v(j-1,i) + v(j  ,i)) * (v(j-1,i) - v(j  ,i))
    )
    ) / (dy*4.0);
  }

  @Inline
  protected def FD_x_UV(
  u: Array[Double]{rank==2},
  v: Array[Double]{rank==2},
  j: Int,
  i: Int,
  dx: Double,
  alpha: Double): Double
  {
    return (           ( (u(j,  i) + u(j+1,  i)) * (v(j,  i) + v(j,i+1)) -
    (u(j,i-1) + u(j+1,i-1)) * (v(j,i-1) + v(j,i  ))
    )
    + alpha * ( Math.abs(u(j,  i) + u(j+1,  i)) * (v(j,i  ) - v(j,i+1)) -
    Math.abs(u(j,i-1) + u(j+1,i-1)) * (v(j,i-1) - v(j,i  ))
    )
    ) / (dx*4.0);
  }

  @Inline
  protected def FD_y_UV(
  u: Array[Double]{rank==2},
  v: Array[Double]{rank==2},
  j: Int,
  i: Int,
  dy: Double,
  alpha: Double): Double
  {
    return (           ( (v(j  ,i) + v(j  ,i+1)) * (u(j  ,i) + u(j+1,i)) -
    (v(j-1,i) + v(j-1,i+1)) * (u(j-1,i) + u(j  ,i))
    )
    + alpha * ( Math.abs(v(j  ,i) + v(j  ,i+1)) * (u(j  ,i) - u(j+1,i)) -
    Math.abs(v(j-1,i) + v(j-1,i+1)) * (u(j-1,i) - u(j  ,i))
    )
    ) / (dy*4.0);
  }

  @Inline
  protected def is_point_in_obstacle_cell(
  M: Array[Int]{rank==2},
  j: Int,
  i: Int) : Boolean
  {
    return ((M(j,i) > 0) || (M(j+1,i) > 0) || (M(j+1,i+1) > 0) || (M(j,i+1) > 0)) ? true : false;
  }

  @Inline
  private def _li(gi: Int): Int {
    return gi - this._gimin + 1 as Int;
  }

  @Inline
  private def _lj(gj: Int): Int {
    return gj - this._gjmin + 1 as Int;
  }

  @Inline
  private def _gi(i: Int): Int {
    return i + this._gimin - 1 as Int;
  }

  @Inline
  private def _gj(j: Int): Int {
    return j + this._gjmin - 1 as Int;
  }

}
