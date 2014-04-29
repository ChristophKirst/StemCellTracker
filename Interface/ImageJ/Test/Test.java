/* Testing Java Stuff */

class Test {
   
   
   /**
    * Class constructor.
    */
   public Test() {
   }

   public void testArray(Object object) {
      System.out.println("test");
      
     if (object instanceof double[][][]) {
        System.out.println("double");
      } else if (object instanceof float[][][]) {
         System.out.println("float");
      } else if (object instanceof short[][][]) {
         System.out.println("short");
      } else if (object instanceof int[][][]) {
         System.out.println("int");
      } else if (object instanceof byte[][][]) {
         System.out.println("byte");
      } else {
         System.out.println("other");
      }

  
      double[][][] is = (double[][][]) object;
      
      int s1 = is.length;
      int s2 = is[0].length;
      int s3 = is[0][0].length;
      System.out.println("dim");
      System.out.println(s1);
      System.out.println(s2);   
      System.out.println(s3);
      
      
      double[][] tt = {{3,4}, {2,3}, {6,7}};
      
      double[] xx = new double[6];
      System.arraycopy(tt[0], 0, xx, 0, 2);
      
      System.out.println("tst");
      System.out.println(xx[0]);
      System.out.println(xx[1]);   
      System.out.println(xx);
      
      
      //float[] yy = (float[]) xx;
      
      //System.out.println(yy[0]);
      //System.out.println(yy[1]);   
      
      
           
   }
   
   public Object get() {
      //byte [] x = new byte[6];
      short [] x = new short[7];
      return x;
   }
   
   
}


