package nssim;

public class Parameters{

  //Boundary types
  public static val BOUNDARY_TYPE_INLET =    -10 as Int;
  public static val BOUNDARY_TYPE_OUTLET =   -20 as Int;
  public static val BOUNDARY_TYPE_NOSLIP =   -30 as Int;
  public static val BOUNDARY_TYPE_FREESLIP = -40 as Int;

  //Simulation types
  public static val SERIAL   = -10 as Int;
  public static val PARALLEL = -20 as Int;
  public static val INVASIVE = -30 as Int;

  public static val NUM_OBS = 17 as Int;
  public static val OBS_XMIN = [0.0, 2.0,  9.0, 11.0, 16.0,  0.0,  5.0, 10.0, 11.0, 16.0, 14.0, 15.0, 18.0, 24.0, 23.0, 22.0, 21.0];
  public static val OBS_XMAX = [9.0, 8.0, 11.0, 16.0, 18.0,  5.0, 10.0, 11.0, 16.0, 21.0, 15.0, 18.0, 19.0, 27.0, 24.0, 23.0, 22.0];
  public static val OBS_YMIN = [0.0, 5.0,  0.0,  0.0,  0.0,  9.0, 10.0, 11.0, 12.0, 13.0,  7.0,  6.0,  8.0,  0.0,  3.0,  4.0,  5.0];
  public static val OBS_YMAX = [5.0, 6.0,  4.0,  3.0,  2.0, 15.0, 15.0, 15.0, 15.0, 15.0,  9.0, 10.0, 10.0, 12.0, 11.0, 10.0, 10.0];

  //Simulation variables
  protected var _domain_size_x: Double = 27;		/// Domain size in x-direction
  protected var _domain_size_y: Double = 15;		/// Domain size in y-direction
  protected var _initial_velocity_x: Double = 1;	/// Initial velocity in x-direction
  protected var _initial_velocity_y: Double = 0;	/// Initial velocity in y-direction
  protected var _initial_pressure: Double = 0;	/// Initial pressure
  protected var _inlet_velocity_x: Double = 1;	/// Inlet velocity in x-direction
  protected var _inlet_velocity_y: Double = 0;	/// Inlet velocity in y-direction
  protected var _external_force_x: Double = 0;	/// External force in x-direction
  protected var _external_force_y: Double = 0;	/// External force in y-direction
  protected var _re: Double = 100;		/// Reynolds number
  protected var _tau: Double = 0.5;		/// Safety factor for time step size computation
  protected var _alpha: Double = 0.9;	/// Upwind differecing factor
  protected var _omega: Double = 1.0;	/// Pressure related
  protected var _boundary_north: Int = BOUNDARY_TYPE_OUTLET;	/// North boundary type
  protected var _boundary_south: Int = BOUNDARY_TYPE_OUTLET;	/// South boundary type
  protected var _boundary_east: Int = BOUNDARY_TYPE_INLET;	/// East boundary type
  protected var _boundary_west: Int = BOUNDARY_TYPE_OUTLET;	/// West boundary type

  //Simulation domain resolution
  protected var _ncx: Int = 27;	/// Local number of cells in x-direction
  protected var _ncy: Int = 15;	/// Local number of cells in y-direction
  protected var _dx: Double;	/// Cell size in x-dimension (same for all blocks)
  protected var _dy: Double;	/// Cell size in y-dimension (same for all blocks)

  // _resx and _resy sets the multiple of base resolution
  // E.g.: _resx = 4 means 27x4 cells in x-direction
  //       _resy = 4 means 15x4 cells in y-direction
  protected var _resx:Int = 2;
  protected var _resy:Int = 2;

  public var final_time:Double = 3.0;
  public var simulation_type:Int = INVASIVE; //SERIAL || PARALLEL || INVASIVE
  public var NP:Int = 4; // Number of processors in the parallel simulation

  // configure the outputs
  public static val vtk_flag:Int = 1;
  public static val output_name:String = "out/artery";


  //MASK VALUES [CEWNS]
  public static val FLUID = 0     as Int;
  public static val B_E   = 10111 as Int;
  public static val B_W   = 11011 as Int;
  public static val B_N   = 11101 as Int;
  public static val B_S   = 11110 as Int;
  public static val B_NE  = 10101 as Int;
  public static val B_SE  = 10110 as Int;
  public static val B_NW  = 11001 as Int;
  public static val B_SW  = 11010 as Int;
  public static val B_IN  = 11111 as Int;

}
