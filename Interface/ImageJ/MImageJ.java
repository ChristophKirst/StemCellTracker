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
    * Create a 2D ImagePlus object from double array
    *
    * @param title    title of the new image  
    * @param object   Matlab image array representing a 2D image
    * @return the resulting ImagePlus instance
    */
   public static ImagePlus createImage2D(String title, Object object) {

         if (object instanceof double[][]) {
            //if (verbose)
            //         System.out.println("MIJ warning message: Loss of precision: convert double 32-bit to float 32-bit");
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
        
         /// Color Images
         else if (object instanceof double[][][]) {
                double[]
                int height = is.length;
                int width = is[0].length;
                int stackSize = is[0][0].length;
                ColorProcessor colorprocessor = new ColorProcessor(width, height);
                byte[] R_pixels = new byte[width * height];
                byte[] G_pixels = new byte[width * height];
                byte[] B_pixels = new byte[width * height];
                if (stackSize >= 3) {
                        for (int h = 0; h < height; h++) {
                                int index = h * width;
                                int w = 0;
                                while (w < width) {
                                        R_pixels[index] = (byte) is[h][w][0];
                                        G_pixels[index] = (byte) is[h][w][1];
                                        B_pixels[index] = (byte) is[h][w][2];
                                        w++;
                                        index++;
                                }
                        }
                } 
                else if (stackSize >= 2) {
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
                else if (stackSize >= 1) {
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
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         else if (object instanceof byte[][][]) {
               byte[][][] is = (byte[][][]) object;
               int height = is.length;
               int width = is[0].length;
               int stackSize = is[0][0].length;
               ImageStack imagestack = new ImageStack(width, height);
               for (int sz = 0; sz < stackSize; sz++) {
                        ByteProcessor byteprocessor = new ByteProcessor(width, height);
                        byte[] bp = (byte[]) byteprocessor.getPixels();
                        i = 0;
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
                        i = 0;
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
               if (verbose)
                        System.out.println("MIJ warning message: Loss of precision: convert int 32 bits to short 16 bits");
               int[][][] is = (int[][][]) object;
               int height = is.length;
               int width = is[0].length;
               int stackSize = is[0][0].length;
               ImageStack imagestack = new ImageStack(width, height);
               for (int sz = 0; sz < stackSize; sz++) {
                        ShortProcessor shortprocessor  = new ShortProcessor(width, height);
                        short[] sp = (short[]) shortprocessor.getPixels();
                        i = 0;
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
                        i = 0;
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
               if (verbose)
                        System.out.println("MIJ warning message: Loss of precision: convert double 32-bit to float 32-bit");
               double[][][] ds = (double[][][]) object;
               int height = ds.length;
               int width = ds[0].length;
               int stackSize = ds[0][0].length;
               ImageStack imagestack = new ImageStack(width, height);
               for (int sz = 0; sz < stackSize; sz++) {
                        FloatProcessor floatprocessor = new FloatProcessor(width, height);
                        float[] fp = (float[]) floatprocessor.getPixels();
                        i = 0;
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
               System.out.println("MIJ Error message: Unknow type of images or volumes.");
               return null;
         }

         if (showImage) {
               imp.show();
               imp.updateAndDraw();
         }
         return imp;
      }

        /**
         * Create a new color image in ImageJ from a Matlab variable.
         *
         * The last index of the array is the color channel index (3 channels) in
         * the follwing order Red-Green-Blue.
         *
         * @param is    Matlab variable 
         */
        public static ImagePlus createColor(final byte[][][] is, boolean showImage) {
                return createColor("Imported from Matlab", is, showImage);
        }

        /**
         * Create a new color image in ImageJ from a Matlab variable with a specified title.
         *
         * @param title title of the new image
         * @param is    Matlab variable
         */
        public static ImagePlus createColor(String title, final byte[][][] is, boolean showImage) {
                int height = is.length;
                int width = is[0].length;
                int stackSize = is[0][0].length;
                ColorProcessor colorprocessor = new ColorProcessor(width, height);
                byte[] R_pixels = new byte[width * height];
                byte[] G_pixels = new byte[width * height];
                byte[] B_pixels = new byte[width * height];
                if (stackSize >= 3) {
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
                else if (stackSize >= 2) {
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
                else if (stackSize >= 1) {
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
                if (showImage) {
                        imp.show();
                }
                return imp;
        }

        /**
         * Create a new 3D color image in ImageJ from a Matlab variable.
         *
         * @param is    Matlab variable
         */
        public static ImagePlus createColor(byte[][][][] is, boolean showImage) {
                return createColor("Import from Matlab", is, showImage);
        }

        /**
         * Create a new 3D color image in ImageJ from a Matlab variable with a specified title.
         *
         * @param title title of the new image
         * @param is    Matlab variable
         */
        public static ImagePlus createColor(String title, byte[][][][] is, boolean showImage) {
                int height = is.length;
                int width = is[0].length;
                int stackSize = is[0][0].length;
                int dim = is[0][0][0].length;
                ImageStack imagestack = new ImageStack(width, height);
                if (dim >= 3) {
                        for (int k = 0; k < stackSize; k++) {
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
                else if (dim >= 2) {
                        for (int k = 0; k < stackSize; k++) {
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
                else if (dim >= 1) {
                        for (int k = 0; k < stackSize; k++) {
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
                if (showImage) {
                        imp.show();
                }
                return imp;
        }
   
   
}
