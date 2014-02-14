/**
 * mij.java
 * Matlab interface to ImageJ and ImagePlus classes
 *
 * Description:
 * Purpose of this class is to provide a basic interface to ImageJ
 * and its ImagePlus class for visualization of 3d images
 * 
 * Chirstoph Kirst
 * The Rockerfeller University, New York, USA
 */

import ij.ImageJ;
import ij.ImagePlus;

class MImageJ {

	public static ImageJ imagej;
   
	/**
	 * Class constructor.
	 */
	public MImageJ() {
      start();
	}

  	/**
	 * start ImageJ
	 */ 
   public void start() {
      if (imagej instanceof ImageJ) return;
      
		imagej = new ImageJ();
    
      if (imagej == null) {
         System.out.println("MImageJ: failed to open ImageJ!");
      }
   }

	/**
	 * Exits ImageJ.imagej instance
	 */
	public static void exit() {
		if (imagej != null) {
         imagej.quit();
         imagej = null;
      }
      
		if (verbose && !(imagej instanceof ImageJ)) {
			System.out.println("ImageJ exited!");
		}
	}

   /**
	 * prints some information about ImageJ settings 
	 */
   public static void info() {
      System.out.println("ImageJ Version:" + IJ.getVersion());
      System.out.println("ImageJ Memory:" + IJ.freeMemory() );
      System.out.println("ImageJ Directory plugins: " + (IJ.getDirectory("plugins")==null?"Not specified":IJ.getDirectory("plugins")));
      System.out.println("ImageJ Directory macros: " + (IJ.getDirectory("macros")==null?"Not specified":IJ.getDirectory("macros")));
      System.out.println("ImageJ Directory luts: " + (IJ.getDirectory("luts")==null?"Not specified":IJ.getDirectory("luts")));
      System.out.println("ImageJ Directory image: " + (IJ.getDirectory("image")==null?"Not specified":IJ.getDirectory("image")));
      System.out.println("ImageJ Directory imagej: " + (IJ.getDirectory("imagej")==null?"Not specified":IJ.getDirectory("imagej")));
      System.out.println("ImageJ Directory startup: " + (IJ.getDirectory("startup")==null?"Not specified":IJ.getDirectory("startup")));
      System.out.println("ImageJ Directory home: " + (IJ.getDirectory("home")==null?"Not specified":IJ.getDirectory("home")));
   }   
}
