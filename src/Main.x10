import nssim.Parameters;
import nssim.NSBlock;
import nssim.NSSim;
import nssim.NSSimInvasive;
import x10.io.File;
import x10.io.FileReader;

public class Main
{
	public static def main(Array[String])
	{
		val param = new Parameters();

		if (param.simulation_type==param.SERIAL) {
			val sim = new NSSim();
			sim.run();
		} else {
			val sim = new NSSimInvasive();
			sim.run();
		}

	}
}
