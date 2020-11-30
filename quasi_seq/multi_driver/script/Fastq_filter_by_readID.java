import java.util.*;
import java.io.*;
import java.util.zip.GZIPInputStream;

public class Fastq_filter_by_readID
{

        public Fastq_filter_by_readID(){}

        //public extract_dogTag(String inputFastq_fileName, String smrtTag, String dogTag_suffix){}

        //public print_statistics(String dogTag_Fastq_fileName){}


    public static void main(String[] args)
    {
    	final int DEFAULT_BUFFER_SIZE=5096;

    	/*String[] names = { "@NCC-1701-D",
           "@NCC-1701-A" };

        Set<String> nameset=new HashSet<String>();
        nameset.addAll(Arrays.asList(names));
*/
        BufferedReader  inputFastq_file=null;

		//BufferedReader inputFastq_file=null;
		BufferedReader inputID_file=null;
		String read_header = null;

        //PrintWriter filtered_fastq = new PrintWriter(new FileOutputStream(args[0]+".readLength_LessThan"+args[1]+".fasta"));
     //   PrintWriter filtered_fastq = new PrintWriter(new FileOutputStream(file_name+".readLength_LessThan"+args[1]+".fasta"));

		String [] input_file_name = args[0].split("/");
        String file_name = input_file_name[0];
        if (input_file_name.length > 0)
        	file_name = input_file_name[input_file_name.length-1];

        System.out.println("input file _name:"+file_name);

		try{

			if(args[0].endsWith(".gz"))
	         {
				 inputFastq_file=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(args[0]),DEFAULT_BUFFER_SIZE)),DEFAULT_BUFFER_SIZE);
				 file_name = file_name.substring(1, file_name.length()-3);
	         }
	     else
	         {
	    	 	inputFastq_file=new BufferedReader(new FileReader(args[0]),DEFAULT_BUFFER_SIZE);
	         }

		String [] idFileName = args[1].split("/");
        String outFileName = idFileName[0];
        if (idFileName.length > 0)
        	outFileName = idFileName[idFileName.length-1];
        PrintWriter filtered_fastq = new PrintWriter(new FileOutputStream(outFileName + ".fastq"));
		//PrintWriter filtered_fastq = new PrintWriter(new FileOutputStream("filtered_"+file_name));
					   // inputFastq_file = new BufferedReader(new FileReader(new File(args[0])));

					    inputID_file = new BufferedReader(new FileReader(new File(args[1])));
					    Hashtable<String, String> readID_hash = new Hashtable<String, String>();

					    while (null != (read_header = inputID_file.readLine()))
					    {
					    	String [] DD = read_header.split(" ");
					    	String [] ID = DD[0].split("/");
					    	if (ID.length>2)
					    		for (int ct=1;ct<ID.length-1;ct++)
					    			ID[0] = ID[0] + "/"+ID[ct];
					    	readID_hash.put("@"+ID[0], "1");
						}
					    while (null != (read_header = inputFastq_file.readLine()))
					    {
					    	String [] ID = read_header.split(" ");
					          if (null != readID_hash.get(ID[0])) // value found in Hash
					          {
					              filtered_fastq.println(read_header);
					              filtered_fastq.println(inputFastq_file.readLine());
					              filtered_fastq.println(inputFastq_file.readLine());
					              filtered_fastq.println(inputFastq_file.readLine());
					           }else{
					        	   inputFastq_file.readLine();
					        	   inputFastq_file.readLine();
					        	   inputFastq_file.readLine();
					           }
				}
					    inputFastq_file.close();
					    inputID_file.close();
					    filtered_fastq.close();
        		}catch (Exception e) {e.printStackTrace();}

    }
}


/*
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.zip.GZIPInputStream;

public class App {
public static void main(String[] args) throws Exception
    {
    final int DEFAULT_BUFFER_SIZE=5096;
    String[] names = { "@NCC-1701-D",
       "@NCC-1701-A" };

    Set<String> nameset=new HashSet<String>();
    nameset.addAll(Arrays.asList(names));

    BufferedReader r;

if(args[0].endsWith(".gz"))
         {
         r=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(args[0]),DEFAULT_BUFFER_SIZE)),DEFAULT_BUFFER_SIZE);
         }
     else
         {
         r=new BufferedReader(new FileReader(args[0]),DEFAULT_BUFFER_SIZE);
         }
    long nLine=-1L;
    String line;
    while((line=r.readLine())!=null)
        {
        nLine++;
        if(nLine%4!=0 || !nameset.contains(line)) continue;
        System.out.println(line);
        System.out.println(r.readLine());
        System.out.println(r.readLine());
        System.out.println(r.readLine());
        nLine+=3;
        }
    r.close();
    }
}

*/
