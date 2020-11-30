/* fileName: ReadCount_snp.java
** Description: use varScan readcount file to generate SNVs
** input: varScan readCount file
** output: SNVs with p-value and frequency
** author: xiaoli.jiao@nih.gov
*/
import java.util.*;
import java.io.*;

public class ReadCount_snp
{
    public ReadCount_snp(){}

    public static void main(String[] args)
    {
    	if (args.length < 1) {
			System.out.println("Usage: java ReadCount_snp xxxx.pileup.readcounts output_filename  estimate_error_rate[default 0.01] p-Value_threshold[default 0.05] varFreq_threshold[default 0.001] coverage_threshold[default 10]");
		    System.exit(0);
    	}
    	double estimate_error_rate = 0.01;
    	double pValue_threshold = 0.05;
    	int coverage_threshold = 10;
	    double varFreq_threshold = 0.001;

    	if (args.length >2)
    		estimate_error_rate = Double.parseDouble(args[2]);
    	if (args.length >3)
    		pValue_threshold = Double.parseDouble(args[3]);
    	if (args.length >4)
    		coverage_threshold = Integer.parseInt(args[4]);
		if (args.length >5)
    		varFreq_threshold = Double.parseDouble(args[5]);


    	final int DEFAULT_BUFFER_SIZE=5096;
        BufferedReader  inputCns_file=null;

        String [] input_file_name = args[0].split("/");
        String file_name = input_file_name[0];

        if (input_file_name.length > 0)
        	file_name = input_file_name[input_file_name.length-1];

        String [] name_str = file_name.split("\\.");
        //String output_file = name_str[0]+".readcounts.snp";

        String output_file = file_name +".snp";
        if (args.length >1)
             output_file = args[1];


		try{
			inputCns_file=new BufferedReader(new FileReader(args[0]),DEFAULT_BUFFER_SIZE);
			PrintWriter snpFile = new PrintWriter(new FileOutputStream(output_file));

			String line = null;
			String [] header = null;
			String fa = "";
			String code = null;
			if (null != (line = inputCns_file.readLine()))	//readcounts header
				header = line.split("\t");

			int pos=0;
			int currPos=0;
			int startPos=0;
			String ins_pos="";
			String insertion="";
			snpFile.println("Chrom\tPosition\tRef\tVarAllele\tReads1\tReads2\tVarFreq\tPvalue");
			while (null != (line = inputCns_file.readLine()))
		    {
				//System.out.print(line+"\n");
				String [] field = line.split("\t");
				if (Integer.parseInt(field[4]) > coverage_threshold)
				{
					pos = Integer.parseInt(field[1]);
					//System.out.println(pos+"  pos \n"+line+"\n");
					if (field.length > 5)
					{
						int depth20=Integer.parseInt(field[4]);
						int cc=0;
						int[] freq = new int[4];
						//String insertion = "";
						int max_ins = 0;
						for (int j=5;j<field.length;j++)
						{
							//System.out.print("field["+j+"]="+field[j]+"\n");
							String [] str = field[j].split(":");
							if (str[0].length() < 2)
							{
								switch(str[0].charAt(0)) {
									case 'A': freq[0]=Integer.parseInt(str[1]);
										if (j==5 && str.length > 7)
										   if (freq[0] < Integer.parseInt(str[7]))
											  freq[0] = Integer.parseInt(str[7]);
									break;
									case 'C': freq[1]=Integer.parseInt(str[1]);
										if (j==5 && str.length > 7 && freq[1] < Integer.parseInt(str[7]))
											freq[1] = Integer.parseInt(str[7]);
									break;
									case 'G': freq[2]=Integer.parseInt(str[1]);
										if (j==5 && str.length > 7 &&  freq[2] < Integer.parseInt(str[7]))
											freq[2] = Integer.parseInt(str[7]);
									break;
									case 'T': freq[3]=Integer.parseInt(str[1]);
										if (j==5 && str.length > 7 &&  freq[3] < Integer.parseInt(str[7]))
											freq[3] = Integer.parseInt(str[7]);
									break;
									default: break;
								}
						/*	}else{  //check if there is an insertion
								if (str[0].startsWith("INS"))
								{
									int read_count=Integer.parseInt(str[1]);

									if (read_count > 7 && read_count > depth20*0.05)
									{	if (cc < 1 )
										{
											//System.out.println();
											for (int i=0; i<5;i++)
												{}//System.out.print(field[i]+"\t");
											//System.out.print(depth20+"\t");
										}
										float pct = read_count*100/depth20;
										//System.out.print(field[j]+"\t"+ pct+"%\t");
										cc=cc+1;
									}
									if (Integer.parseInt(str[1])> max_ins )
									{
										max_ins = Integer.parseInt(str[1]);
										String [] ins = str[0].split("-");  //INS-6-ACATGG
										insertion = ins[2]; //ACATGG
									}
								}

								if (str[0].startsWith("DEL"))
								{
									int read_count=Integer.parseInt(str[1]);
									if (read_count > 7 && read_count > depth20*0.05)
									{	if (cc < 1 )
										{
											System.out.println();
											for (int i=0; i<5;i++)
												System.out.print(field[i]+"\t");
										}
										float pct = read_count*100/depth20;
										//System.out.print(field[j]+"\t"+ pct+"%\t");
										cc=cc+1;
									}
								}*/
							}
						}

						int freq_sum = freq[0]+freq[1]+freq[2]+freq[3];

						if (freq_sum < coverage_threshold)
						    continue;

						if (freq_sum>depth20*0.3)
						{	//continue;

							int max_count=0;int second_max_count=0;int max_idx=-1; int second_max_idx=-1;

							for (int i=0;i<4;i++)
								if (max_count < freq[i])
								{
									max_count = freq[i];
									max_idx = i;
								}
							for (int i=0;i<4;i++)
								if (i == max_idx)
									continue;
								else
									if (second_max_count < freq[i])
									{
										second_max_count =freq[i];
										second_max_idx = i;
									}

						//float varFreq = new float;
							double varFreq = second_max_count*100.0/freq_sum;
							String refAllele="N";
							String varAllele="N";
							if (varFreq > varFreq_threshold*100)
							{
								//System.out.println("\nmax_idx="+max_idx+" second_max_idx="+second_max_idx);

								switch (max_idx) {
								//case 0: return "-"+insertion;
									case 0: refAllele="A";break;
									case 1: refAllele="C";break;
									case 2: refAllele="G";break;
									case 3: refAllele="T";break;
									default: break;
								}
								switch (second_max_idx) {
								//case 0: return "-"+insertion;
									case 0: varAllele="A";break;
									case 1: varAllele="C";break;
									case 2: varAllele="G";break;
									case 3: varAllele="T";break;
									default: break;
								}
								FishersExact fe = new FishersExact(100000);
								double estimate_errors = max_count*estimate_error_rate;
								//System.out.print("estimate errors:\t"+ estimate_errors);
								double pValue = fe.getRightTailedP(max_count, (int)estimate_errors, max_count, second_max_count);
								if (pValue < pValue_threshold)
								   snpFile.println(field[0]+"\t"+field[1]+"\t"+refAllele+"\t"+varAllele+"\t"+max_count+"\t"+second_max_count+"\t"+varFreq+"\t"+pValue);
							}
						}
					}
				}
		    }

			inputCns_file.close();
			//inputRef_fasta.close();
			snpFile.close();

		}catch (Exception e) {e.printStackTrace();}

    }

}
