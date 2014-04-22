/**
 * MImageJ.java
 * Matlab interface to ImageJ and ImagePlus classes
 *
 * Description:
 * Purpose of this class is to provide a basic interface to ImageJ
 * and its ImagePlus class for visualization of 3d images
 * 
 * Christoph Kirst
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
      System.out.println("ImageJ Free Memory:" + IJ.freeMemory() );
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
    * @param foramt   image format a string like 'pq', 'pql' 'pqlc' 'pqltc' (p, q, l = pixel coordinates, c = color, t= time)
    * @return the ImagePlus
    */
   public static ImagePlus createImage(String title, String format, Object object) {
      ImagePlus img;
      if (format.equals("pq")) {
         img = createImageGray2D(title, object);
      } else if (format.equals("pql")) {
         img = createImageGray3D(title, object);
      } else if (format.equals("lpq")) {
         img = createImageGray3DFast(title, object);              
      } else if (format.equals("pqc")) {
         img = createImageRGB2D(title, object);
      } else if (format.equals("cpq")) {
         img = createImageRGB2DFast(title, object);
      } else if (format.equals("pqlc")) {
         img = createImageRGB3D(title, object);
      } else if (format.equals("clpq")) {
         img = createImageRGB3DFast(title, object); 
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
			ByteProcessor byteprocessor = new ByteProcessor(width, height);
			byte[] bp = (byte[]) byteprocessor.getPixels();
			for (int h = 0; h < height; ++h) {
            System.arraycopy(is[h], 0, bp, h*width, width);
         }
			imp = new ImagePlus(title, byteprocessor);
		} 
      
		else if (object instanceof short[][]) {
			short[][] is = (short[][]) object;
			int height = is.length;
			int width = is[0].length;
			ShortProcessor shortprocessor = new ShortProcessor(width, height);
			short[] sp = (short[]) shortprocessor.getPixels();
			for (int h = 0; h < height; ++h) {
            System.arraycopy(is[h], 0, sp, h*width, width);
         }
			imp = new ImagePlus(title, shortprocessor);
		} 
      
		else if (object instanceof int[][]) {
			int[][] is = (int[][]) object;
			int height = is.length;
			int width = is[0].length;
			ShortProcessor shortprocessor = new ShortProcessor(width, height);
			short[] sp = (short[]) shortprocessor.getPixels();
			for (int h = 0; h < height; ++h) {
            System.arraycopy(is[h], 0, sp, h*width, width);
         }
			imp = new ImagePlus(title, shortprocessor);
		} 
      
		else if (object instanceof float[][]) { 
			float[][] fs = (float[][]) object;
			int height = fs.length;
			int width = fs[0].length;
			FloatProcessor floatprocessor = new FloatProcessor(width, height);
			float[] fp = (float[])floatprocessor.getPixels();
			for (int h = 0; h < height; ++h) {
            System.arraycopy(fs[h], 0, fp, h*width, width);
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
			for (int h = 0; h < height; ++h) {
            System.arraycopy(ds[h], 0, fp, h*width, width);
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
    * @param object   Matlab image array representing a gray 3D image (pql format)
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
				for (int h= 0; h < height; ++h) {
					for (int w = 0; w < width; w++) {
						bp[i] = is[h][w][sz];
						i++;
					}
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
				for (int h= 0; h < height; ++h) {
					for (int w = 0; w < width; w++) {
						sp[i] = is[h][w][sz];
						i++;
					}
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
				for (int h= 0; h < height; ++h) {
					for (int w = 0; w < width; w++) {
						sp[i] = (short) is[h][w][sz];
						i++;
					}
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
				for (int h = 0; h < height; ++h) {
					for (int w = 0; w < width; ++w) {
						fp[i] = fs[h][w][sz];
						i++;
					}
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
				for (int h= 0; h < height; ++h) {
					for (int w = 0; w < width; w++) {
						fp[i] = (float) ds[h][w][sz];
						i++;
					}
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
    * Create ImagePlus object from a gray 3d Matlab array (fast)
    *
    * @param title    title of the new image  
    * @param object   Matlab image array representing a gray 3D image (lpq format)
    * @return the ImagePlus
    */
   public static ImagePlus createImageGray3DFast(String title, Object object) {
      // switch types

      ImagePlus imp;
      
		if (object instanceof byte[][][]) {
			byte[][][] is = (byte[][][]) object;
			int stackSize = is.length;
			int height = is[0].length;
			int width = is[0][0].length;
			ImageStack imagestack = new ImageStack(width, height);
			for (int sz = 0; sz < stackSize; sz++) {
				ByteProcessor byteprocessor = new ByteProcessor(width, height);
				byte[] bp = (byte[]) byteprocessor.getPixels();
            for (int h = 0; h < height; ++h) {
               System.arraycopy(is[sz][h], 0, bp, h*width, width);
            }
				imagestack.addSlice("", byteprocessor);
			}
			imp = new ImagePlus(title, imagestack);
		} 
      
		else if (object instanceof short[][][]) {
			short[][][] is = (short[][][]) object;
			int stackSize = is.length;
			int height = is[0].length;
			int width = is[0][0].length;
			ImageStack imagestack = new ImageStack(width, height);
			for (int sz = 0; sz < stackSize; sz++) {
				ShortProcessor shortprocessor  = new ShortProcessor(width, height);
				short[] sp = (short[]) shortprocessor.getPixels();
            for (int h = 0; h < height; ++h) {
               System.arraycopy(is[sz][h], 0, sp, h*width, width);
            }
				imagestack.addSlice("", shortprocessor);
			}
			imp = new ImagePlus(title, imagestack);
		} 
      
		else if (object instanceof int[][][]) {
			int[][][] is = (int[][][]) object;
			int stackSize = is.length;
			int height = is[0].length;
			int width = is[0][0].length;
			ImageStack imagestack = new ImageStack(width, height);
			for (int sz = 0; sz < stackSize; sz++) {
				ShortProcessor shortprocessor  = new ShortProcessor(width, height);
				short[] sp = (short[]) shortprocessor.getPixels();
				int i = 0;
				for (int h= 0; h < height; ++h) {
					for (int w = 0; w < width; w++) {
						sp[i] = (short) is[h][w][sz];
						i++;
					}
				}
				if (sz == 0)
					shortprocessor.resetMinAndMax();
				imagestack.addSlice("", shortprocessor);
			}
			imp = new ImagePlus(title, imagestack);
		} 
      
		else if (object instanceof float[][][]) {
			float[][][] fs = (float[][][]) object;
			int stackSize = fs.length;
			int height = fs[0].length;
			int width = fs[0][0].length;
			ImageStack imagestack = new ImageStack(width, height);
			for (int sz = 0; sz < stackSize; sz++) {
				FloatProcessor floatprocessor = new FloatProcessor(width, height);
				float[] fp = (float[]) floatprocessor.getPixels();
            for (int h = 0; h < height; ++h) {
               System.arraycopy(fs[sz][h], 0, fp, h*width, width);
            }
				if (sz == 0) floatprocessor.resetMinAndMax();
				imagestack.addSlice("", floatprocessor);
			}
			imp=new ImagePlus(title, imagestack);
		} 
      
		else if (object instanceof double[][][]) {
			double[][][] ds = (double[][][]) object;
			int stackSize = ds.length;
			int height = ds[0].length;
			int width = ds[0][0].length;
			ImageStack imagestack = new ImageStack(width, height);
			for (int sz = 0; sz < stackSize; sz++) {
				FloatProcessor floatprocessor = new FloatProcessor(width, height);
				float[] fp = (float[]) floatprocessor.getPixels();
				int i = 0;
				for (int h= 0; h < height; ++h) {
					for (int w = 0; w < width; w++) {
						fp[i] = (float) ds[h][w][sz];
						i++;
					}
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
    * Create ImagePlus object from a color 2D Matlab array (pqc)
    *
    * @param title    title of the new image  
    * @param object   Matlab image array representing a rgb 2D image (pqc format)
    * @return the ImagePlus
    */         
    public static ImagePlus createImageRGB2D(String title, Object object) {
       
      if (!(object instanceof byte[][][])) {
         System.out.println("MImageI: error: image is not a byte (unit8) RGB 2D array.");
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
         int index = 0;
			for (int h = 0; h < height; h++) {
            for (int w = 0; w < width; w++) {
					R_pixels[index] = is[h][w][0];
					G_pixels[index] = is[h][w][1];
					B_pixels[index] = is[h][w][2];
					index++;
				}
			}
		} 
		else if (colsize == 2) {
         int index = 0;
			for (int h = 0; h < height; h++) {
            for (int w = 0; w < width; w++) {
					R_pixels[index] = is[h][w][0];
					G_pixels[index] = is[h][w][1];
					index++;
				}
			}
		} 
		else if (colsize == 1) {
         int index = 0;
			for (int h = 0; h < height; h++) {
            for (int w = 0; w < width; w++) {
					R_pixels[index] = is[h][w][0];
					index++;
				}
			}
		}
      
		colorprocessor.setRGB(R_pixels, G_pixels, B_pixels);
		ImagePlus imp = new ImagePlus(title, colorprocessor);

		return imp;
	}

    
   /**
    * Create ImagePlus object from a color 2D Matlab array (fast, cpq)
    *
    * @param title    title of the new image  
    * @param object   Matlab image array representing a rgb 2D image (cpq format)
    * @return the ImagePlus
    */         
    public static ImagePlus createImageRGB2DFast(String title, Object object) {
       
      if (!(object instanceof byte[][][])) {
         System.out.println("MImageI: error: image is not a byte (unit8) RGB 2D array.");
			return null;
		}

      byte[][][] is = (byte[][][]) object;

		int colsize = is.length;
		int height = is[0].length;
		int width = is[0][0].length;
		ColorProcessor colorprocessor = new ColorProcessor(width, height);
		byte[] R_pixels = new byte[width * height];
		byte[] G_pixels = new byte[width * height];
		byte[] B_pixels = new byte[width * height];
		if (colsize >= 3) {
			for (int h = 0; h < height; h++) {
            System.arraycopy(is[0][h], 0, R_pixels, h * width, width);
            System.arraycopy(is[1][h], 0, G_pixels, h * width, width);  
            System.arraycopy(is[2][h], 0, B_pixels, h * width, width);
			}
		} 
		else if (colsize == 2) {
			for (int h = 0; h < height; h++) {
            System.arraycopy(is[0][h], 0, R_pixels, h * width, width);
            System.arraycopy(is[1][h], 0, G_pixels, h * width, width);  
			}
		} 
		else if (colsize == 1) {
			for (int h = 0; h < height; h++) {
            System.arraycopy(is[0][h], 0, R_pixels, h * width, width);
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
    * @param object   Matlab image array representing a rgb 3D image (pqlc format)
    * @return the ImagePlus
    */     
	public static ImagePlus createImageRGB3D(String title, Object object) {
      
      if (!(object instanceof byte[][][][])) {
         System.out.println("MImageI: error: image is not a byte (unit8) RGB 3D array.");
			return null;
		}
      
      byte[][][][] is = (byte[][][][]) object;
      
		int height = is.length;
		int width = is[0].length;
		int length = is[0][0].length;
		int colsize = is[0][0][0].length;
		ImageStack imagestack = new ImageStack(width, height);
		if (colsize >= 3) {
			for (int l = 0; l < length; l++) {
				ColorProcessor colorprocessor = new ColorProcessor(width, height);
				byte[] R_pixels = new byte[width * height];
				byte[] G_pixels = new byte[width * height];
				byte[] B_pixels = new byte[width * height];   
            int index = 0;
				for (int h = 0; h < height; h++) {
               for (int w = 0; w < width; w++) {
						R_pixels[index] = is[h][w][l][0];
						G_pixels[index] = is[h][w][l][1];
						B_pixels[index] = is[h][w][l][2];
						index++;
					}
				}
				colorprocessor.setRGB(R_pixels, G_pixels, B_pixels);
				imagestack.addSlice("", colorprocessor);
			}
		} 
		else if (colsize == 2) {
			for (int l = 0; l < length; l++) {
				ColorProcessor colorprocessor = new ColorProcessor(width, height);
				byte[] R_pixels = new byte[width * height];
				byte[] G_pixels = new byte[width * height];
				byte[] B_pixels = new byte[width * height];   
            int index = 0;
				for (int h = 0; h < height; h++) {
               for (int w = 0; w < width; w++) {
						R_pixels[index] = is[h][w][l][0];
						G_pixels[index] = is[h][w][l][1];
						index++;
					}
				}
				colorprocessor.setRGB(R_pixels, G_pixels, B_pixels);
				imagestack.addSlice("", colorprocessor);
			}
		} 
		else if (colsize == 1) {
			for (int l = 0; l < length; l++) {
				ColorProcessor colorprocessor = new ColorProcessor(width, height);
				byte[] R_pixels = new byte[width * height];
				byte[] G_pixels = new byte[width * height];
				byte[] B_pixels = new byte[width * height];   
            int index = 0;
				for (int h = 0; h < height; h++) {
               for (int w = 0; w < width; w++) {
						R_pixels[index] = is[h][w][l][0];
						index++;
					}
				}
				colorprocessor.setRGB(R_pixels, G_pixels, B_pixels);
				imagestack.addSlice("", colorprocessor);
			}
		}
		ImagePlus imp = new ImagePlus(title, imagestack);

		return imp;
	}
   
   
   
   /**
    * Create ImagePlus object from a rgb 3d Matlab array (fast, clpq format)
    *
    * @param title    title of the new image  
    * @param object   Matlab image array representing a rgb 3D image (clpq format)
    * @return the ImagePlus
    */     
	public static ImagePlus createImageRGB3DFast(String title, Object object) {
      
      if (!(object instanceof byte[][][][])) {
         System.out.println("MImageI: error: image is not a byte (unit8) RGB 3D array.");
			return null;
		}
      
      byte[][][][] is = (byte[][][][]) object;
      
		int colsize = is.length;
		int length = is[0].length;
		int height = is[0][0].length;
		int width = is[0][0][0].length;
      
      //System.out.println(colsize);
      //System.out.println(length);
      //System.out.println(height);
      //System.out.println(width);
      
		ImageStack imagestack = new ImageStack(width, height);
		if (colsize >= 3) {
			for (int l = 0; l < length; l++) {
				ColorProcessor colorprocessor = new ColorProcessor(width, height);
				byte[] R_pixels = new byte[width * height];
				byte[] G_pixels = new byte[width * height];
				byte[] B_pixels = new byte[width * height];   
            int index = 0;
				for (int h = 0; h < height; h++) {
               System.arraycopy(is[0][l][h], 0, R_pixels, h * width, width);
               System.arraycopy(is[1][l][h], 0, G_pixels, h * width, width);
               System.arraycopy(is[2][l][h], 0, B_pixels, h * width, width);
				}
				colorprocessor.setRGB(R_pixels, G_pixels, B_pixels);
				imagestack.addSlice("", colorprocessor);
			}
		} 
		else if (colsize == 2) {
			for (int l = 0; l < length; l++) {
				ColorProcessor colorprocessor = new ColorProcessor(width, height);
				byte[] R_pixels = new byte[width * height];
				byte[] G_pixels = new byte[width * height];
				byte[] B_pixels = new byte[width * height];   
            int index = 0;
				for (int h = 0; h < height; h++) {
               System.arraycopy(is[0][l][h], 0, R_pixels, h * width, width);
               System.arraycopy(is[1][l][h], 0, G_pixels, h * width, width);
				}
				colorprocessor.setRGB(R_pixels, G_pixels, B_pixels);
				imagestack.addSlice("", colorprocessor);
         }
		} 
		else if (colsize == 1) {
			for (int l = 0; l < length; l++) {
				ColorProcessor colorprocessor = new ColorProcessor(width, height);
				byte[] R_pixels = new byte[width * height];
				byte[] G_pixels = new byte[width * height];
				byte[] B_pixels = new byte[width * height];   
            int index = 0;
				for (int h = 0; h < height; h++) {
               System.arraycopy(is[0][l][h], 0, R_pixels, h * width, width);
				}
				colorprocessor.setRGB(R_pixels, G_pixels, B_pixels);
				imagestack.addSlice("", colorprocessor);
         }
		}
		ImagePlus imp = new ImagePlus(title, imagestack);

		return imp;
	}
   
 
} // class MImageJ
