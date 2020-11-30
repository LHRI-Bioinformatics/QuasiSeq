//getClust.java
import com.mathworks.toolbox.javabuilder.*;
import SigClust.*;

class getClust
{
   public static void main(String[] args)
   {
    //  MWNumericArray n = null;
    //  Object[] result = null;
      Class1 theMagic = null;

      if (args.length < 3)
      {
        System.out.println("Error: must input 3 strings: 1. xxx.bam 2.snv_pos.txt 3.prefix of output headers");
        return;
      }

      try
      {
        // n = new MWNumericArray(Double.valueOf(args[0]),
               //                       MWClassID.DOUBLE);
          String s1 = new String(args[0]);
          String s2 = new String(args[1]);
          String s3 = new String(args[2]);

          //String s1 = new String("Hello");
          //String s2 = new String("there");
          //String s3 = new String("Vandana");

          MWCharArray x = new MWCharArray(s1);
          MWCharArray y = new MWCharArray(s2);
          MWCharArray z = new MWCharArray(s3);

          //System.out.println(x);
          //System.out.println(y);
          //System.out.println(z);
         theMagic = new Class1();
		theMagic.SigClust(x,y,z);
        // result = theMagic.zpClust(1, n);
        // System.out.println(result[0]);
      }
      catch (Exception e)
      {
         System.out.println("Exception: " + e.toString());
      }
      finally
      {
        // MWArray.disposeArray(n);
        // MWArray.disposeArray(result);
	if (theMagic != null)
         theMagic.dispose();
      }
   }
}
