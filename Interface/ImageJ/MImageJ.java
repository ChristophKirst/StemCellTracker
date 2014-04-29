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
   
   
   
   /**
    * Returns the current (selected) image from ImageJ.
    *
    * @return	Current image
    */
   public static Object getCurrentImage() {
      ImagePlus imageplus = WindowManager.getCurrentImage();
      return fromImagePlus(imageplus);
   }
   
   /**
    * Returns the specified image from ImageJ.
    *
    * @param title	title of image
    * @return Object
    */
   public static Object getImage(String title) {
      String[] strings = getListImages();
      int[] is = WindowManager.getIDList();
      for (int i = 0; i < is.length; i++) {
         if (strings[i].equals(title)) {
            ImagePlus imageplus = WindowManager.getImage(is[i]);
            return get(imageplus);
         }
      }
      System.out.println("MIJ Error message: the requested image (" + title + ") does not exist.");
      return null;
   }


   
   
   /**
    * Convert ImagePlus to Array
    *
    * @param imageplus	image
    * @return an N x M array representing the input image
    */
   private static Object fromImagePlus(ImagePlus imageplus) {
      if (imageplus == null)
         return null;
      int width = imageplus.getWidth();
      int height = imageplus.getHeight();
      int length = imageplus.getStackSize();
      ImageStack imagestack = imageplus.getStack();

      switch (imageplus.getType()) {

         case ImagePlus.GRAY8: {
            byte[][][] is = new byte[height][width][stackSize];
            for (int l = 0; l < length; l++) {
               ByteProcessor byteprocessor = (ByteProcessor) imagestack.getProcessor(l + 1);
               byte[] pixels = (byte[]) byteprocessor.getPixels();
               int counter = 0;
               for (int h = 0: h < height; h++) {
                  for (int w = 0; w < width; w++) {
                     is[h][w][l] = pixels[counter];
                     counter++;
                  }
               }
            }
            return is;
         }
         case ImagePlus.GRAY16: {
            short[][][] is = new short[height][width][length];
            for (int l = 0; l < length; l++) {
               ShortProcessor shortprocessor = (ShortProcessor) imagestack.getProcessor(sz + 1);
               short[] pixels = (short[]) shortprocessor.getPixels();
               int counter = 0;
               for (int h = 0: h < height; h++) {
                  for (int w = 0; w < width; w++) {
                     is[h][w][l] = pixels[counter];
                     counter++;
                  }
               }
            }
            return is;
         }
         case ImagePlus.GRAY32: {
            int[][][] fs = new int[height][width][stackSize];
            for (int sz = 0; sz < stackSize; sz++) {
               FloatProcessor floatprocessor = (FloatProcessor) imagestack.getProcessor(sz + 1);
               float[] fpixels = (float[]) floatprocessor.getPixels();
               counter = 0;
               int i = 0;
               while (i < height) {
                  int j = 0;
                  while (j < width) {
                     fs[i][j][sz] = (double) fpixels[counter];
                     j++;
                     counter++;
                  }
                  counter = ++i * width;
               }
            }
            return fs;
         }
         
         
         
         
                  
         case ImagePlus.COLOR_256: {
            ;
         }
         
         
         
         case ImagePlus.COLOR_RGB: {
            if (stackSize == 1) {
               short[][][] is = new short[height][width][3];
               ColorProcessor colorprocessor = (ColorProcessor) imagestack.getProcessor(1);
               byte[] red = new byte[width * height];
               byte[] green = new byte[width * height];
               byte[] blue = new byte[width * height];
               colorprocessor.getRGB(red, green, blue);
               counter = 0;
               int h = 0;
               while (h < height) {
                  int w = 0;
                  while (w < width) {
                     is[h][w][0] = (short)(red[counter]&0xff);
                     is[h][w][1] = (short)(green[counter]&0xff);
                     is[h][w][2] = (short)(blue[counter]&0xff);
                     w++;
                     counter++;
                  }
                  counter = ++h * width;
               }
               return is;
            }
            short[][][][] is = new short[height][width][stackSize][3];
            for (int sz = 0; sz < stackSize; sz++) {
               ColorProcessor colorprocessor  = (ColorProcessor) imagestack.getProcessor(sz + 1);
               byte[] red = new byte[width * height];
               byte[] green = new byte[width * height];
               byte[] blue = new byte[width * height];
               colorprocessor.getRGB(red, green, blue);
               counter = 0;
               int h = 0;
               while (h < height) {
                  int w = 0;
                  while (w < width) {
                     is[h][w][sz][0] = (short)red[counter];
                     is[h][w][sz][1] = (short)green[counter];
                     is[h][w][sz][2] = (short)blue[counter];
                     w++;
                     counter++;
                  }
                  counter = ++h * width;
               }
            }
            return is;
         }
         default:
            System.out.println("MIJ Error message: Unknow type of volumes.");
            return null;
      }
   }
   
   
   
   
   
   
   
   
   
   
   
   /**
    * Gives the list of the open images in ImageJ.
    *
    * @return	List of open images names
    */
   public static String[] getImageTitleList() {
      int[] is = WindowManager.getIDList();
      String[] strings = new String[is.length];
      for (int i = 0; i < is.length; i++) {
         ImagePlus imageplus = WindowManager.getImage(is[i]);
         strings[i] = imageplus.getTitle();
      }
      return strings;
   }
   
   /**
    * Returns the title of the current image window.
    *
    * @return	Title of the current image window
    */
   public static String getCurrentImageTitle() {
      ImagePlus imageplus = WindowManager.getCurrentImage();
      return imageplus.getTitle();
   }
   

   
   
   
   

   
   
   
   
   
   
   
   
   
   
   
   
   
   /**
    * Set a region of interest (ROI) in the current image.
    *
    * @param	roiarray	give coordinates or positions of the ROI depending of the ROI type
    * @param	type	supported types: Roi.LINE, Roi.RECTANGLE, Roi.POINT, Roi.OVAL, Roi.POLYLINE, Roi.POLYGON, Roi.ANGLE
    */
   public static void setRoi(double[][] roiarray, int type) {
      ImagePlus imageplus = WindowManager.getCurrentImage();
      switch (type) {
         case Roi.LINE:
            Line roi=new Line((int) roiarray[0][0], (int) roiarray[1][0],
            (int) roiarray[0][1], (int) roiarray[1][1]);
            imageplus.setRoi((Roi)roi);
            break;
         case Roi.RECTANGLE:
            int width= (int) roiarray[0][0]-(int) roiarray[0][1];
            int height= (int) roiarray[1][1]-(int) roiarray[1][2];
            Roi rect=new Roi((int) roiarray[0][0], (int) roiarray[1][0], Math.abs(width), Math.abs(height));
            imageplus.setRoi(rect);
            break;
         case Roi.POINT:
            PointRoi pnt=new PointRoi((int) roiarray[0][0], (int) roiarray[1][0], imageplus);
            imageplus.setRoi(pnt);
            break;
         case Roi.OVAL:
            width= (int) roiarray[0][0]-(int) roiarray[0][1];
            height= (int) roiarray[1][1]-(int) roiarray[1][2];
            OvalRoi oval=new OvalRoi((int) roiarray[0][0], (int) roiarray[1][0], Math.abs(width), Math.abs(height));
            imageplus.setRoi(oval);
            break;
         case Roi.POLYLINE:
            int nPoints=roiarray[0].length;
            int[] xarr=new int[nPoints];
            int[] yarr=new int[nPoints];
            for (int i=0;i<nPoints;i++) {
               xarr[i]=(int)roiarray[0][i];
               yarr[i]=(int)roiarray[1][i];
            }
            PolygonRoi poly=new PolygonRoi(xarr, yarr, nPoints, Roi.POLYLINE);
            imageplus.setRoi(poly);
            break;
         case Roi.POLYGON:
            nPoints=roiarray[0].length;
            xarr=new int[nPoints];
            yarr=new int[nPoints];
            for (int i=0;i<nPoints;i++) {
               xarr[i]=(int)roiarray[0][i];
               yarr[i]=(int)roiarray[1][i];
            }
            poly=new PolygonRoi(xarr, yarr, nPoints, Roi.POLYGON);
            imageplus.setRoi(poly);
            break;
         case Roi.ANGLE:
            break;
         default:
      }
   }
   
   /**
    * Get a region of interest (ROI) of the current image with or without calibration.
    *
    * @param option	CAL for using calibration or NOCAL for no calibration
    * @return Object
    */
   public static Object getRoi(int option) {
      ImagePlus imageplus = WindowManager.getCurrentImage();
      Roi roi=imageplus.getRoi();
      Calibration cal=imageplus.getCalibration();
      double fh=cal.pixelHeight;
      double fw=cal.pixelWidth;
      Object ret=null;
      if (roi.isLine()) {
         Rectangle rect=roi.getBounds();
         double x=rect.getX();
         double y=rect.getY();
         double w=rect.getWidth();
         double h=rect.getHeight();
         if (option==NOCAL){
            double [][] pnts= {{x,x+w,x+w,x},{y,y,y+h,y+h}};
            ret=(Object)pnts;
         }
         if (option==CAL){
            double [][] pnts= {{x*fw,(x+w)*fw,(x+w),x*fw},
            {y*fh,y*fh,(y+h)*fh,(y+h)*fh}};
            ret=(Object)pnts;
         }
      }
      else {
         Polygon  polygon=roi.getPolygon();
         if (option==NOCAL){
            int[][] pnts=new int[2][polygon.npoints];
            pnts[0]=polygon.xpoints;
            pnts[1]=polygon.ypoints;
            ret=(Object)pnts;
         }
         if (option==CAL){
            double [][] pnts=new double[2][polygon.npoints];
            for (int i=0;i<polygon.npoints; i++){
               pnts[0][i]=polygon.xpoints[i]*fw;
               pnts[1][i]=polygon.ypoints[i]*fh;
            }
            ret=(Object)pnts;
         }
      }
      return ret;
   }
   

   
   
   
   
   
      
   /**
    * Exits ImageJ.imagej instance
    */
   public static void exit() {
      IJ.getInstance().quit();
      imagej = null;
      if (imagej instanceof ImageJ) {
         System.out.println("MImageJ: could not quit ImageJ instance!");
      }
   }
   
 
} // class MImageJ
