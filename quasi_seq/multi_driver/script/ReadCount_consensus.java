/* fileName: ReadCount_consensus
** Description: use varScan readcount file to generate the consensus sequence with ambiguity codes
** input: varScan readCount file
** output: consensus sequence with ambiguityCode
** author: xiaoli.jiao@nih.gov
*/
import java.util.*;
import java.io.*;

public class ReadCount_consensus
{
    public ReadCount_consensus(){}

    /*	ambiguityCode:
 	M          A or C
    R          A or G
    W          A or T
    S          C or G
    Y          C or T
    K          G or T
    V        A or C or G
    H        A or C or T
    D        A or G or T
    B        C or G or T
    X      G or A or T or C
    */
    public static String get_ambiguityCode(String [] field, int ambiguity_threshold, int coverage)
    {
		if (field.length > 5)
		{
			int depth20=Integer.parseInt(field[4]);
			int cc=0;
			int[] freq = new int[4];
			String insertion = "";
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
				}else{  //check if there is an insertion
					if (str[0].startsWith("INS"))
					{
						int read_count=Integer.parseInt(str[1]);

						if (read_count > 7 && read_count > depth20*0.05)
						{	if (cc < 1 )
							{
								System.out.println();
								for (int i=0; i<5;i++)
									System.out.print(field[i]+"\t");
								//System.out.print(depth20+"\t");
							}
							float pct = read_count*100/depth20;
							System.out.print(field[j]+"\t"+ pct+"%\t");
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
							System.out.print(field[j]+"\t"+ pct+"%\t");
							cc=cc+1;
						}
					}
				}
			}

			int freq_sum = freq[0]+freq[1]+freq[2]+freq[3];
			if (freq_sum==0)
			    freq_sum = depth20;

			int max_count=0;
			for (int j=0;j<4;j++)
				if (max_count < freq[j])
					max_count = freq[j];
			if	(depth20 < coverage && max_count < 8)
			{	return "N";
			}else{	//if (max_count > 7)
				float[] percent = new float[4];
				String aboveThr = "";
				for (int i = 0; i < 4; i++)
				{
					percent[i] = freq[i]*100/freq_sum;
					if (percent[i] >= ambiguity_threshold)
						aboveThr = aboveThr + "1";
					else
						aboveThr = aboveThr + "0";
				}

				if (max_ins < Integer.parseInt(field[4])*0.5)
					insertion = "";

				if (null != aboveThr)
				{
					switch (Integer.parseInt(aboveThr)) {
					//case 0: return "-"+insertion;
					case 0: return insertion;
					case 1000: return "A"+insertion;
					case 100: return "C"+insertion;
					case 10: return "G"+insertion;
					case 1: return "T"+insertion;
					case 1100: return "M"+insertion;
					case 1010: return "R"+insertion;
					case 1001: return "W"+insertion;
					case 110: return "S"+insertion;
					case 101: return "Y"+insertion;
					case 11: return "K"+insertion;
					case 1110: return "V"+insertion;
					case 1101: return "H"+insertion;
					case 1011: return "D"+insertion;
					case 111: return "B"+insertion;
					case 1111: return "X"+insertion;
					default: return null;
					}
				}else {return null;}
			}
		}
		return null;
    }

    public static void main(String[] args)
    {
    	if (args.length < 1) {
			System.out.println("Usage: java ReadCount_consensus xxxx.pileup.readcounts ambiguity_threshold[default 50] coverage_threshold[default 8]");
		    System.exit(0);
    	}
    	int abiguity_threshold = 50;
    	int coverage_threshold = 8;
    	if (args.length >1)
    		abiguity_threshold = Integer.parseInt(args[1]);
    	if (args.length >2)
    		coverage_threshold = Integer.parseInt(args[2]);

    	final int DEFAULT_BUFFER_SIZE=5096;
        BufferedReader  inputCns_file=null;

        String [] input_file_name = args[0].split("/");
        String file_name = input_file_name[0];

        if (input_file_name.length > 0)
        	file_name = input_file_name[input_file_name.length-1];

        String [] name_str = file_name.split("\\.");


        String output_file = name_str[0]+".consensus.fasta";

		try{
			inputCns_file=new BufferedReader(new FileReader(args[0]),DEFAULT_BUFFER_SIZE);
			PrintWriter cns_fasta = new PrintWriter(new FileOutputStream(output_file));

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

			while (null != (line = inputCns_file.readLine()))
		    {
				System.out.print(line+"\n");
				String [] field = line.split("\t");
				if (Integer.parseInt(field[4]) > coverage_threshold){
					pos = Integer.parseInt(field[1]);
					startPos = pos;
					if (null != (code = get_ambiguityCode(field,abiguity_threshold,coverage_threshold)))
					    fa = code;
					else
						fa = "";

					break;
				}
		    }

			while (null != (line = inputCns_file.readLine()))
		    {
				String [] field = line.split("\t");
				currPos = Integer.parseInt(field[1]);
				if (currPos - pos > 1 )
				{
					for (int k=pos+1; k<currPos;k++)
						fa = fa +"N";

					if (null != (code = get_ambiguityCode(field,abiguity_threshold,coverage_threshold)))
					{
						if (code.length()>2)
						{
							ins_pos = ins_pos + ","+field[1];
							insertion = insertion+","+code;
						}
						fa = fa + code;
					}

				}else{ //(currPos - pos == 1 )
					if (null != (code = get_ambiguityCode(field,abiguity_threshold,coverage_threshold)))
					{
						if (code.length()>2)
						{
							ins_pos = ins_pos + ","+field[1];
							insertion = insertion+","+code;
						}
						fa = fa + code;
				   }
				}
				pos = currPos;
		     }

			if (pos > 0)
			{
				while (fa.endsWith("N")) //trim end"NNNNNNNNNNNNN"
					fa = fa.substring(0,fa.length()-1);
				System.out.println("\n>"+name_str[0] + ".consensus" + "\tlength="+fa.length()+"\tref:"+startPos+"-"+pos + "\tInsertions at:"+ins_pos + "\t"+insertion);
				cns_fasta.println(">"+name_str[0]+".consensus"+ " length="+fa.length());
				cns_fasta.println(fa);
			}else{
				System.out.println("\n>"+name_str[0] + "\tinsufficient coverage");
				cns_fasta.println(">"+name_str[0]+".consensus"  + " length=0 \t insufficient coverage" );

			}

			inputCns_file.close();
			//inputRef_fasta.close();
			cns_fasta.close();

		}catch (Exception e) {e.printStackTrace();}

    }

}
