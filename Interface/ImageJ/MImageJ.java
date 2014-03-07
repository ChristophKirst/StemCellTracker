/**
 * MImageJ.java
 * Matlab interface to ImageJ and ImagePlus classes
 *
 * Description:
 * Purpose of this class is to provide a basic interface to ImageJ
 * and its ImagePlus class for visualization of 3d images
 * 
 * Chirstoph Kirst
 * The Rockerfeller University, New York, USA
 *
 * based on MIJ.java
 */


import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ShortProcessor;


class MImageJ {

   /**
    * Class constructor.
    */
   public MImageJ() {
   }

   /**
    * Prints information about ImageJ settings 
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
   
   /**
    * Create ImagePlus object for Matlab array using a psecified format
    *
    * @param title    title of the new image  
    * @param object   Matlab image array representing a 2D image
    * @param foramt   image format a string like 'hw', 'hwl' 'hwlc' 'hwltc' (h = height, w = width, l = length, c = color, t= time)
    * @return the ImagePlus
    */
   public static ImagePlus createImage(String title, String format, Object object) {
      ImagePlus img;
      if (format.equals("hw")) {
         img = createImageGray2D(title, object);
      } else if (format.equals("hwl")) {
         img = createImageGray3D(title, object);
      } else if (format.equals("hwc")) {
         img = createImageRGB2D(title, object);
      } else if (format.equals("hwlc")) {
         img = createImageRGB3D(title, object);
      } else{
         System.out.println("MImageI: error: unknown image format: " + format);
         return null;
      }
      
      return img;
   }

   /**
    * Create ImagePlus object from a gray 2d Matlab array 
    *
    * @param title    title of the new image  
    * @param object   Matlab image array representing a 2D gray image (hw format)
    * @return the ImagePlus
    */
   public static ImagePlus createImageGray2D(String title, Object object) {
      //switch between types
      int i = 0;
      ImagePlus imp;
      
		if (object instanceof byte[][]) {
			byte[][] is = (byte[][]) object;
			int height = is.length;
			int width = is[0].length;
			ByteProcessor byteprocessor = new ByteProcessor(width, height); // check ordering here !! width height ??
			byte[] bp = (byte[]) byteprocessor.getPixels();
			int h = 0;
			while (h < height) {
				int w = 0;
				while (w < width) {
					bp[i] = is[h][w];
					w++;
					i++;
				}
				i = ++h * width;
			}
			imp = new ImagePlus(title, byteprocessor);

		} 
      
		else if (object instanceof short[][]) {
         
			short[][] is = (short[][]) object;
			int height = is.length;
			int width = is[0].length;
			ShortProcessor shortprocessor = new ShortProcessor(width, height);
			short[] sp = (short[]) shortprocessor.getPixels();
			int h = 0;
			while (h < height) {
				int w = 0;
				while (w < width) {
					sp[i] = is[h][w];
					w++;
					i++;
				}
				i = ++h * width;
			}
			imp = new ImagePlus(title, shortprocessor);

		} 
      
		else if (object instanceof int[][]) {

			int[][] is = (int[][]) object;
			int height = is.length;
			int width = is[0].length;
			ShortProcessor shortprocessor = new ShortProcessor(width, height);
			short[] sp = (short[]) shortprocessor.getPixels();
			int h = 0;
			while (h < height) {
				int w = 0;
				while (w < width) {
					sp[i] = (short)is[h][w];
					w++;
					i++;
				}
				i = ++h * width;
			}
			imp = new ImagePlus(title, shortprocessor);
		} 
      
		else if (object instanceof float[][]) {
         
			float[][] fs = (float[][]) object;
			int height = fs.length;
			int width = fs[0].length;
			FloatProcessor floatprocessor = new FloatProcessor(width, height);
			float[] fp = (float[])floatprocessor.getPixels();
			int h = 0;
			while (h < height) {
				int w = 0;
				while (w < width) {
					fp[i] = fs[h][w];
					w++;
					i++;
				}
				i = ++h * width;
			}
			floatprocessor.resetMinAndMax();
			imp = new ImagePlus(title, floatprocessor);

		} 
      
		else if (object instanceof double[][]) {
         
			double[][] ds = (double[][]) object;
			int height = ds.length;
			int width = ds[0].length;
			FloatProcessor floatprocessor = new FloatProcessor(width, height);
			float[] fp = (float[]) floatprocessor.getPixels();
			int h = 0;
			while (h < height) {
				int w = 0;
				while (w < width) {
					fp[i] = (float) ds[h][w];
					w++;
					i++;
				}
				i = ++h * width;
			}
			floatprocessor.resetMinAndMax();
			imp = new ImagePlus(title, floatprocessor);

		} 
      
		else {
			System.out.println("MImageI: error: image is not of format BW2D.");
			return null;
		}

		return imp;
	}
           
  
   /**
    * Create ImagePlus object from a gray 3d Matlab array 
    *
    * @param title    title of the new image  
    * @param object   Matlab image array representing a gray 3D image (hwl format)
    * @return the ImagePlus
    */
   public static ImagePlus createImageGray3D(String title, Object object) {
      // switch types

      ImagePlus imp;
      
		if (object instanceof byte[][][]) {
			byte[][][] is = (byte[][][]) object;
			int height = is.length;
			int width = is[0].length;
			int stackSize = is[0][0].length;
			ImageStack imagestack = new ImageStack(width, height);
			for (int sz = 0; sz < stackSize; sz++) {
				ByteProcessor byteprocessor = new ByteProcessor(width, height);
				byte[] bp = (byte[]) byteprocessor.getPixels();
				int i = 0;
				int h = 0;
				while (h < height) {
					int w = 0;
					while (w < width) {
						bp[i] = is[h][w][sz];
						w++;
						i++;
					}
					i = ++h * width;
				}
				imagestack.addSlice("", byteprocessor);
			}
			imp = new ImagePlus(title, imagestack);

		} 
      
		else if (object instanceof short[][][]) {
			short[][][] is = (short[][][]) object;
			int height = is.length;
			int width = is[0].length;
			int stackSize = is[0][0].length;
			ImageStack imagestack = new ImageStack(width, height);
			for (int sz = 0; sz < stackSize; sz++) {
				ShortProcessor shortprocessor  = new ShortProcessor(width, height);
				short[] sp = (short[]) shortprocessor.getPixels();
				int i = 0;
				int h = 0;
				while (h < height) {
					int w = 0;
					while (w < width) {
						sp[i] = is[h][w][sz];
						w++;
						i++;
					}
					i = ++h * width;
				}
				imagestack.addSlice("", shortprocessor);
			}
			imp = new ImagePlus(title, imagestack);

		} 
      
		else if (object instanceof int[][][]) {
			int[][][] is = (int[][][]) object;
			int height = is.length;
			int width = is[0].length;
			int stackSize = is[0][0].length;
			ImageStack imagestack = new ImageStack(width, height);
			for (int sz = 0; sz < stackSize; sz++) {
				ShortProcessor shortprocessor  = new ShortProcessor(width, height);
				short[] sp = (short[]) shortprocessor.getPixels();
				int i = 0;
				int h = 0;
				while (h < height) {
					int w = 0;
					while (w < width) {
						sp[i] = (short) is[h][w][sz];
						w++;
						i++;
					}
					i = ++h * width;
				}
				if (sz == 0)
					shortprocessor.resetMinAndMax();
				imagestack.addSlice("", shortprocessor);

			}
			imp = new ImagePlus(title, imagestack);

		} 
      
		else if (object instanceof float[][][]) {
			float[][][] fs = (float[][][]) object;
			int height = fs.length;
			int width = fs[0].length;
			int stackSize = fs[0][0].length;
			ImageStack imagestack = new ImageStack(width, height);
			for (int sz = 0; sz < stackSize; sz++) {
				FloatProcessor floatprocessor = new FloatProcessor(width, height);
				float[] fp = (float[]) floatprocessor.getPixels();
				int i = 0;
				int h = 0;
				while (h < height) {
					int w = 0;
					while (w < width) {
						fp[i] = fs[h][w][sz];
						w++;
						i++;
					}
					i = ++h * width;
				}
				if (sz == 0)
					floatprocessor.resetMinAndMax();
				imagestack.addSlice("", floatprocessor);
			}
			imp=new ImagePlus(title, imagestack);

		} 
      
		else if (object instanceof double[][][]) {
			double[][][] ds = (double[][][]) object;
			int height = ds.length;
			int width = ds[0].length;
			int stackSize = ds[0][0].length;
			ImageStack imagestack = new ImageStack(width, height);
			for (int sz = 0; sz < stackSize; sz++) {
				FloatProcessor floatprocessor = new FloatProcessor(width, height);
				float[] fp = (float[]) floatprocessor.getPixels();
				int i = 0;
				int h = 0;
				while (h < height) {
					int w = 0;
					while (w < width) {
						fp[i] = (float) ds[h][w][sz];
						w++;
						i++;
					}
					i = ++h * width;
				}
				if (sz == 0)
					floatprocessor.resetMinAndMax();
				imagestack.addSlice("", floatprocessor);
			}
			imp=new ImagePlus(title, imagestack);

		}
      
      else {
			System.out.println("MImageI: error: image is not of format BW3D.");
			return null;
		}

      return imp;
   }


   /**
    * Create ImagePlus object from a gray 3d Matlab array 
    *
    * @param title    title of the new image  
    * @param object   Matlab image array representing a rgb 2D image (hwc format)
    * @return the ImagePlus
    */         
    public static ImagePlus createImageRGB2D(String title, Object object) {
       
      if (! ((object instanceof double[][][]) || (object instanceof float[][][]) ||
            (object instanceof short[][][]) || (object instanceof int[][][]) ||
            (object instanceof byte[][][]))) {
         System.out.println("MImageI: error: image is not of format RGB2D.");
			return null;
		}

      byte[][][] is = (byte[][][]) object;

		int height = is.length;
		int width = is[0].length;
		int colsize = is[0][0].length;
		ColorProcessor colorprocessor = new ColorProcessor(width, height);
		byte[] R_pixels = new byte[width * height];
		byte[] G_pixels = new byte[width * height];
		byte[] B_pixels = new byte[width * height];
		if (colsize >= 3) {
			for (int h = 0; h < height; h++) {
				int index = h * width;
				int w = 0;
				while (w < width) {
					R_pixels[index] = is[h][w][0];
					G_pixels[index] = is[h][w][1];
					B_pixels[index] = is[h][w][2];
					w++;
					index++;
				}
			}
		} 
		else if (colsize >= 2) {
			for (int j = 0; j < height; j++) {
				int index = j * width;
				int i = 0;
				while (i < width) {
					R_pixels[index] = is[j][i][0];
					G_pixels[index] = is[j][i][1];
					i++;
					index++;
				}
			}
		} 
		else if (colsize >= 1) {
			for (int j = 0; j < height; j++) {
				int index = j * width;
				int i = 0;
				while (i < width) {
					R_pixels[index] = is[j][i][0];
					i++;
					index++;
				}
			}
		}
      
		colorprocessor.setRGB(R_pixels, G_pixels, B_pixels);
		ImagePlus imp = new ImagePlus(title, colorprocessor);

		return imp;
	}

   /**
    * Create ImagePlus object from a rgb 3d Matlab array 
    *
    * @param title    title of the new image  
    * @param object   Matlab image array representing a rgb 3D image (hwlc format)
    * @return the ImagePlus
    */     
	public static ImagePlus createImageRGB3D(String title, Object object) {
      
       if (! ((object instanceof double[][][][]) || (object instanceof float[][][][]) ||
            (object instanceof short[][][][]) || (object instanceof int[][][][]) ||
            (object instanceof byte[][][][]))) {
         System.out.println("MImageI: error: image is not of format RGB3D.");
			return null;
		}

      byte[][][][] is = (byte[][][][]) object;
      
		int height = is.length;
		int width = is[0].length;
		int length = is[0][0].length;
		int colsize = is[0][0][0].length;
		ImageStack imagestack = new ImageStack(width, height);
		if (colsize >= 3) {
			for (int k = 0; k < length; k++) {
				ColorProcessor colorprocessor = new ColorProcessor(width, height);
				byte[] red = new byte[width * height];
				byte[] green = new byte[width * height];
				byte[] blue = new byte[width * height];
				for (int j = 0; j < height; j++) {
					int index = j * width;
					int i = 0;
					while (i < width) {
						red[index] = is[j][i][k][0];
						green[index] = is[j][i][k][1];
						blue[index] = is[j][i][k][2];
						i++;
						index++;
					}
				}
				colorprocessor.setRGB(red, green, blue);
				imagestack.addSlice("", colorprocessor);
			}
		} 
		else if (colsize >= 2) {
			for (int k = 0; k < length; k++) {
				ColorProcessor colorprocessor = new ColorProcessor(width, height);
				byte[] red = new byte[width * height];
				byte[] green = new byte[width * height];
				byte[] blue = new byte[width * height];

				for (int j = 0; j < height; j++) {
					int index = j * width;
					int i = 0;
					while (i < width) {
						red[index] = is[j][i][k][0];
						green[index] = is[j][i][k][1];
						i++;
						index++;
					}
				}
				colorprocessor.setRGB(red, green, blue);
				imagestack.addSlice("", colorprocessor);
			}
		} 
		else if (colsize >= 1) {
			for (int k = 0; k < length; k++) {
				ColorProcessor colorprocessor = new ColorProcessor(width, height);
				byte[] red = new byte[width * height];
				byte[] green = new byte[width * height];
				byte[] blue = new byte[width * height];

				for (int j = 0; j < height; j++) {
					int index = j * width;
					int i = 0;
					while (i < width) {
						red[index] = is[j][i][k][0];
						i++;
						index++;
					}
				}
				colorprocessor.setRGB(red, green, blue);
				imagestack.addSlice("", colorprocessor);
			}
		}
		ImagePlus imp = new ImagePlus(title, imagestack);

		return imp;
	}
}
