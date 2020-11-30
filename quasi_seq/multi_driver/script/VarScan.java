/**
 * @(#)VarScan.java
 *
 * Copyright (c) 2009-2013 Daniel C. Koboldt and Washington University in St. Louis
 *
 * COPYRIGHT
 */

package net.sf.varscan;

//Import required packages //

import java.io.*;
import java.util.*;
import java.text.*;


/**
 * A set of tools for variant detection in next-generation sequence data.
 *
 * @version	2.4
 *
 * @author Daniel C. Koboldt <dkoboldt@genome.wustl.edu>
 *
 * <BR>
 * <pre>
 * COMMANDS
 * pileup2snp [pileup file] OPTIONS
 * 			Call SNPs from a pileup file that meet certain cutoffs
 * 			Input: 	Pileup file and parameters
 * 			Output: SNPs file with read counts and p-value
 *
 * pileup2indel [pileup file] OPTIONS
 * 			Call indels from a pileup file that meet certain cutoffs
 * 			Input: 	Pileup file and parameters
 * 			Output: Indels file with read counts and p-value
 *
 * pileup2cns [pileup file] OPTIONS
 * 			Call consensus genotypes (reference or variant) at sites with sufficient coverage
 * 			Input: 	Pileup file and parameters
 * 			Output: Consensus file with genotypes, read counts and p-values
 *
 * mpileup2cns [pileup file] OPTIONS
 * 			Call consensus genotypes (reference or variant) across one or more samples
 * 			Input: 	SAMtools mpileup file and parameters
 * 			Output: Consensus file with genotypes, read counts and p-values, or VCF file
 *
 * somatic [normal_pileup] [tumor_pileup] [output] OPTIONS
 * 			Determine somatic status of SNPs from normal/tumor pileup for positions
 * 			Input:	Normal pileup, tumor pileup, and positions file
 * 			Output: SNPs file with read counts and somatic status
 *
 * readcounts [pileup] --variants-file [positions] --output-file [output]
 * 			Obtain read counts for each allele of variants from a pileup file
 * 			Input:	Variants file and pileupfile
 * 			Output: Variants file with read counts for each allele
 *
 * filter [variant file] OPTIONS
 * 			Filter a set of SNPs/indels based on coverage, reads, p-value, etc.
 * 			Input:	SNPs file with read counts and p-value
 * 			Output: Filtered SNPs file with read counts and p-value
 *
 * somaticFilter [somatic-status file] OPTIONS
 * 			Filter VarScan Somatic/Germline/LOH calls for clusters and proximal indels
 * 			Input:	VarScan output for SNPs or Indels (varscan.output.snp)
 * 			Output: Variants passing all filters (varscan.output.snp.filter)
 *
 * fpFilter [variant-file] [readcount-file] OPTIONS
 * 			Apply the false-positive filter to VarScan variant calls
 * 			Input:	VarScan output for SNPs or Indels (varscan.output.snp) and bam-readcount output
 * 			Output: Variants passing all filters (varscan.output.snp.fpfilter)
 *
 * processSomatic [somatic-status file] OPTIONS
 * 			Process VarScan output by somatic status and confidence
 * 			Input:	VarScan output for SNPs or Indels (varscan.output.snp)
 * 			Output: Variants by somatic status (varscan.output.snp.Somatic)
 *
 * copyCaller [copynumber file] OPTIONS
 * 			Process VarScan copynumber output to adjust for GC and make preliminary calls
 * 			Input:	VarScan copynumber output (varscan.output.copynumber)
 * 			Output: Normalized copy number with preliminary calls (varscan.output.copynumber.called)
 *
 * compare [file1] [file2] [type] [output] OPTIONS
 * 			Compares chromosome-position entries in two tab-delimited files
 * 			Input:	File 1 and File 2
 * 			Output: Merged, intersected, or unique entries
 *
 * limit [variants] --regions-file [regions] --output-file [output]
 * 			Limit a tab-delimited file (SNPs, pileup, etc) to a set of positions or regions
 * 			Input:	tab-delimited input file with chromosome & position; positions-file or regions-file
 * 			Output: Entries in input-file matching regions or positions
 *
 * coverage [pileup-file] --regions-file [regions] --output-file [output]
 * 			**Experimental** Calculate Q>20 coverage depth/breadth for a set of target regions
 * 			Input:	Pileup file and tab-delimited regions-file
 * 			Output: Coverage report at various Q>20 depths (1x,10x,20x...)

 *
 * </pre>
 *
 *
 */
public class VarScan {

	final static double MIN_FREQ_FOR_HOM = 0.70;

	/**
	 * Runs the main execution logic
	 * @param args		Command-line arguments
	 */
	public static void main(String[] args) {

		String usage = "VarScan v2.4.2\n\n***NON-COMMERCIAL VERSION***\n\nUSAGE: java -jar VarScan.jar [COMMAND] [OPTIONS] \n\n";
		usage = usage + "COMMANDS:\n" +
				"\tpileup2snp\t\tIdentify SNPs from a pileup file\n" +
				"\tpileup2indel\t\tIdentify indels a pileup file\n" +
				"\tpileup2cns\t\tCall consensus and variants from a pileup file\n" +

				"\tmpileup2snp\t\tIdentify SNPs from an mpileup file\n" +
				"\tmpileup2indel\t\tIdentify indels an mpileup file\n" +
				"\tmpileup2cns\t\tCall consensus and variants from an mpileup file\n\n" +

				"\tsomatic\t\t\tCall germline/somatic variants from tumor-normal pileups\n" +
				"\tcopynumber\t\t\tDetermine relative tumor copy number from tumor-normal pileups\n" +
				"\treadcounts\t\tObtain read counts for a list of variants from a pileup file\n\n" +

				"\tfilter\t\t\tFilter SNPs by coverage, frequency, p-value, etc.\n" +
				"\tsomaticFilter\t\tFilter somatic variants for clusters/indels\n" +
				"\tfpfilter\t\tApply the false-positive filter\n\n" +
				"\tprocessSomatic\t\tIsolate Germline/LOH/Somatic calls from output\n" +
				"\tcopyCaller\t\tGC-adjust and process copy number changes from VarScan copynumber output\n" +

				"\tcompare\t\t\tCompare two lists of positions/variants\n" +
				"\tlimit\t\t\tRestrict pileup/snps/indels to ROI positions\n" +
				"\n";

		if(args.length > 0)
		{
			HashMap<String, String> params = getParams(args);

			if(args[0].equals("pileup2snp"))
			{
				pileup2call(args, params, "SNP");
			}

			else if(args[0].equals("pileup2indel"))
			{
				pileup2call(args, params, "INDEL");
			}

			else if(args[0].equals("pileup2cns"))
			{
				pileup2call(args, params, "CNS");
			}

			else if(args[0].equals("mpileup2snp") || args[0].equals("mpileup2indel") || args[0].equals("mpileup2cns") || args[0].equals("mpileup2vcf"))
			{
				mpileup2call(args, params, "CNS");
			}


			else if(args[0].equals("filter"))
			{
				filter(args, params);
			}

			else if(args[0].equals("somaticFilter"))
			{
				somaticFilter(args, params);
			}

			else if(args[0].equals("fpfilter"))
			{
				fpfilter(args, params);
			}

			else if(args[0].equals("processSomatic"))
			{
				processSomatic(args, params);
			}

			else if(args[0].equals("copyCaller"))
			{
				copyCaller(args, params);
			}

			else if(args[0].equals("compare"))
			{
				compare(args, params);
			}

			else if(args[0].equals("readcounts"))
			{
				readcounts(args, params);
			}

			else if(args[0].equals("somatic"))
			{
				somatic(args, params);
			}

			else if(args[0].equals("trio"))
			{
				trio(args, params, "CNS");
			}

			else if(args[0].equals("copynumber"))
			{
				copynumber(args, params);
			}

			else if(args[0].equals("limit"))
			{
				limit(args, params);
			}
			else if(args[0].equals("coverage"))
			{
				coverage(args, params);
			}
			else if(args[0].equals("test"))
			{
				System.err.println("Testing...");
				try
				{
					RandomAccessFile ref = new RandomAccessFile("test.fasta", "r");
					ref.seek(52);
					byte[] buffer = new byte[5];
					ref.read(buffer);
					String thisBase = new String(buffer);
					System.err.println("Got " + thisBase);
				}
				catch(Exception e)
				{
					System.err.println("Error: Reference file: " + e.getLocalizedMessage());
				}

			}

			else
			{
				System.err.println("Command not recognized\n" + usage);
			}
		}
		else
		{
			System.err.println(usage);
		}

	}


	/**
	 * Calls SNPs from a pileup file
	 *
	 * @param	args			Command-line arguments and parameters
	 * @param	callType		"SNP", "INDEL", or "CNS"
	 */
	public static void pileup2call(String[] args, HashMap<String, String> params, String callType)
	{
		CallPileup pileupCall = new CallPileup(args, callType);
	}

	/**
	 * Calls SNPs from an mpileup file
	 *
	 * @param	args			Command-line arguments and parameters
	 * @param	callType		"SNP", "INDEL", or "CNS"
	 */
	public static void mpileup2call(String[] args, HashMap<String, String> params, String callType)
	{
		CallMpileup mpileupCall = new CallMpileup(args, callType);
	}

	/**
	 * Obtains read counts for a list of variants
	 *
	 * @param	args	Command-line arguments
	 */
	public static void readcounts(String[] args, HashMap<String, String> params)
	{
		 ReadCounts myReadCounts = new ReadCounts(args, params);
	}


	/**
	 * Calls somatic/germline/LOH variants from normal and tumor pileup files
	 *
	 * @param	args	Command-line arguments
	 */
	public static void somatic(String[] args, HashMap<String, String> params)
	{
		if(params.containsKey("mpileup"))
		{
			Somatic mySomatic = new Somatic(args, true);
		}
		else
		{
			Somatic mySomatic = new Somatic(args);
		}

	}

	/**
	 * Calls SNPs in a father-mother-child trio from an mpileup file
	 *
	 * @param	args			Command-line arguments and parameters
	 * @param	callType		"SNP", "INDEL", or "CNS"
	 */
	public static void trio(String[] args, HashMap<String, String> params, String callType)
	{
		Trio myTrio = new Trio(args, callType);
	}


	/**
	 * Determines tumor copy number from normal and tumor pileup files
	 *
	 * @param	args	Command-line arguments
	 */
	public static void copynumber(String[] args, HashMap<String, String> params)
	{
		if(params.containsKey("mpileup"))
		{
			 Copynumber myCopynumber = new Copynumber(args, true);
		}
		else
		{
			 Copynumber myCopynumber = new Copynumber(args);
		}

	}


	/**
	 * Filters variants by coverage, significance, frequency, etc.
	 *
	 * @param	args	Command-line arguments
	 */
	public static void filter(String[] args, HashMap<String, String> params)
	{
		FilterVariants myFilter = new FilterVariants(args);
	}


	/**
	 * Applies false positive filter using bam-readcount information
	 *
	 * @param	args	Command-line arguments
	 */
	public static void fpfilter(String[] args, HashMap<String, String> params)
	{
		FpFilter myFilter = new FpFilter(args);
	}

	/**
	 * Filters variants by coverage, significance, frequency, etc.
	 *
	 * @param	args	Command-line arguments
	 */
	public static void somaticFilter(String[] args, HashMap<String, String> params)
	{
		FilterSomatic myFilter = new FilterSomatic(args);
	}

	/**
	 * Splits VarScan output according to somatic status and confidence
	 *
	 * @param	args	Command-line arguments
	 */
	public static void processSomatic(String[] args, HashMap<String, String> params)
	{
		ProcessSomatic myProcess = new ProcessSomatic(args);
	}

	/**
	 * Calls somatic copy number events from copynumber output
	 *
	 * @param	args	Command-line arguments
	 */
	public static void copyCaller(String[] args, HashMap<String, String> params)
	{
		CopyCaller myCopy = new CopyCaller(args, params);
	}

	/**
	 * Compares two lists of positions/variants
	 *
	 * @param	args	Command-line arguments
	 */
	public static void compare(String[] args, HashMap<String, String> params)
	{
		Comparison myComparison = new Comparison(args);
	}


	/**
	 * Limits pileup or variant files to a list of positions or regions
	 *
	 * @param	args	Command-line arguments
	 */
	public static void limit(String[] args, HashMap<String, String> params)
	{
		LimitVariants myLimit = new LimitVariants(args);
	}

	/**
	 * Reports region coverage from a BAM file
	 *
	 * @param	args	Command-line arguments
	 */
	public static void coverage(String[] args, HashMap<String, String> params)
	{
		Coverage myCoverage = new Coverage(args);
	}



	/**
	 * Parses and verifies any command-line parameters
	 *
	 * @param	args	Command-line arguments
	 * @return			HashMap of parameter names and their values
	 */
	static HashMap getParams(String[] args)
	{
		HashMap<String, String> params = new HashMap<String, String>();

		// Parse out command line arguments //

		String arg = "";
		String value = "";
		int i = 0, j = 0;

		// Go through each argument in the command line //

		while (i < args.length)
		{
			j = i + 1;
			arg = args[i];

			// If the argument starts with a hyphen, make use of it //

			if (arg.startsWith("-"))
			{
				// Remove leading hyphens //
				while(arg.startsWith("-"))
				{
					arg = arg.replaceFirst("-", "");
				}

				// Parse out parameters followed by values //

				if (i < args.length && j < args.length && !args[j].startsWith("-"))
				{
					value = args[j];
					params.put(arg, value);
				}

				// Set other parameters to true //

				else
				{
					params.put(arg, "true");
				}
			}

			i++;
		}

		return(params);
	}


	/**
	 * Gets the infile from command line or input buffer
	 *
	 * @param	args	Command-line arguments
	 * @return			HashMap of parameter names and their values
	 */
	static BufferedReader getInfile(String[] args)
	{
		BufferedReader in = null;

	    try
	    {
	    	// Declare file-parsing variables //

	    	String line;

	    	// Check for file on command line //

	    	if(args.length > 1 && !args[1].startsWith("-"))
	    	{
	    		File infile = new File(args[1]);
	    		if(infile.exists())
	    		{
	    			// Parse the infile //
	    			System.err.println("Reading input from " + args[1]);
	    			in = new BufferedReader(new SmartFileReader(args[1]));
	    		}
	    		else
	    		{
//    				System.err.println("File not found: " + args[1] + "\n");
//    				System.exit(10);
	    		}
	    	}

	    	// If no file from command line was parsed, try for piped input //

	    	if(in == null)
	    	{
		    	// Check the input stream //
		    	InputStreamReader instream = new InputStreamReader(System.in);
		    	Thread.sleep(1000);

		    	int num_naps = 0;

	    		while(!instream.ready())
	    		{
	    			System.err.println("Input stream not ready, waiting for 5 seconds...");
	    			Thread.sleep(5000);
	    			num_naps++;

	    			if(num_naps >= 100)
	    			{
	    				System.err.println("ERROR: Gave up waiting after 500 seconds...\n");
	    				System.exit(10);
	    			}
	    		}

		    	// If we have piped input, proceed with it //

		    	if(instream.ready())
		    	{
		    		System.err.println("Reading input from STDIN");
			    	in = new BufferedReader(instream);
		    	}
	    	}
	    }
	    catch(Exception e)
	    {
	    	System.err.println("ERROR: Unable to open input stream\n");
	    	System.exit(10);
	    }

		return(in);
	}



	/**
	 * Counts the number, quality, and strands of each allele from a pileup
	 *
	 * @param	refBase		Reference base at this position
	 * @param	readBases	String of read bases from pileup
	 * @param	readQuals	String of read base qualities from pileup
	 * @param	minAvgQual	Integer of minimum required base quality to count a base.
	 * @return	results		HashMap<String, String> of results for each allele
	 */
	static HashMap<String, String> getReadCounts(String refBase, String readBases, String readQuals, int minAvgQual, String mapQuals)
	{
		HashMap<String, Integer> readCounts = new HashMap<String, Integer>();
		HashMap<String, Integer> readCountsPlus = new HashMap<String, Integer>();
		HashMap<String, Integer> readCountsMinus = new HashMap<String, Integer>();
		HashMap<String, Integer> qualitySum = new HashMap<String, Integer>();
		HashMap<String, Integer> mapQualitySum = new HashMap<String, Integer>();
		HashMap<String, String> strandsSeen = new HashMap<String, String>();

		int reads1 = 0;
		int reads1indel = 0;
		String readBase = "";
		String prevBase = "";
		String nextBase = "";
		int baseQuality = 0;
		int prevBaseQuality = 0;
		int mapQuality = 1;
		String strand = "";

		String[] arrBases = readBases.split("");
//		char[] arrBases = readBases.toCharArray();
		char[] arrQualities = readQuals.toCharArray();
		char[] mapQualities = mapQuals.toCharArray();

		// Set booleans for read Start //

		boolean readStart = false;

		// Set quality position offset //
		int j = 0;

		// Go through each base //

		for(int i = 0; i < arrBases.length; i++)
		{
			readBase = arrBases[i];
			if(i == 0 && readBase.length() == 0)
			{
				i++;
				readBase = arrBases[i];
			}

			// Record previous and next base //
			prevBase = "";
			if(i > 1 && i < (arrBases.length - 1))
				prevBase = arrBases[i - 1];

			if(j > 1 && j < (arrQualities.length - 1))
				prevBaseQuality = arrQualities[j - 1] - 33;

			nextBase = "";
			if(i < (arrBases.length - 1))
				nextBase = arrBases[i + 1];

			// Get the quality score //
			if(j < arrQualities.length)
				baseQuality = arrQualities[j] - 33;

			// Get the map quality score //
			if(j < mapQualities.length)
				mapQuality = mapQualities[j] - 33;

//			System.err.println("Got " + readBase + " with quality " + arrQualities[j] + "=" + baseQuality + " Next " + nextBase);
			// A period or comma NOT followed by indel represents a reference base //
			if((readBase.equals(".") || readBase.equals(",")) && !(nextBase.equals("-") || nextBase.equals("+")))
			{
				strand = "+";
				if(readBase.equals(","))
					strand = "-";

				if(baseQuality >= minAvgQual)
				{
					reads1++;

					// Count the strands seen //

					if(strandsSeen.containsKey("ref"))
					{
						String alreadySeen = strandsSeen.get("ref");
						if(!(alreadySeen.length() >= 2 || alreadySeen.equals(strand)))
						{
							strandsSeen.put("ref", (strandsSeen.get("ref") + strand));
						}
					}
					else
					{
						strandsSeen.put("ref", strand);
					}

					// Count strand-based read count //
					if(strand.equals("+"))
					{
						// Plus strand //
						if(readCountsPlus.containsKey("ref"))
						{
							readCountsPlus.put("ref", (readCountsPlus.get("ref") + 1));
						}
						else
						{
							readCountsPlus.put("ref", 1);
						}
					}
					else
					{
						// Minus Strand //
						if(readCountsMinus.containsKey("ref"))
						{
							readCountsMinus.put("ref", (readCountsMinus.get("ref") + 1));
						}
						else
						{
							readCountsMinus.put("ref", 1);
						}
					}

					// Add the quality to the sum //

					if(qualitySum.containsKey("ref"))
					{
						qualitySum.put("ref", (qualitySum.get("ref") + baseQuality));
						mapQualitySum.put("ref", (mapQualitySum.get("ref") + mapQuality));
					}
					else
					{
						qualitySum.put("ref", baseQuality);
						mapQualitySum.put("ref", mapQuality);
					}
				}

				j++;

				readStart = false;
			}
			// SNP Processing //
			else if(readBase.toUpperCase().equals("A") || readBase.toUpperCase().equals("C") || readBase.toUpperCase().equals("G") || readBase.toUpperCase().equals("T"))
			{
				strand = "+";

				if(readBase.equals("a") || readBase.equals("c") || readBase.equals("g") || readBase.equals("t"))
					strand = "-";

				readBase = readBase.toUpperCase();

				// Check that we're not at start or end of read //
				if(baseQuality >= minAvgQual)// && !readStart && !nextBase.equals("$"))
				{
					// Count the read //
					if(readCounts.containsKey(readBase))
					{
						readCounts.put(readBase, (readCounts.get(readBase) + 1));
					}
					else
					{
						readCounts.put(readBase, 1);
					}

					// Count strand-based read count //
					if(strand.equals("+"))
					{
						// Plus strand //
						if(readCountsPlus.containsKey(readBase))
						{
							readCountsPlus.put(readBase, (readCountsPlus.get(readBase) + 1));
						}
						else
						{
							readCountsPlus.put(readBase, 1);
						}
					}
					else
					{
						// Minus Strand //
						if(readCountsMinus.containsKey(readBase))
						{
							readCountsMinus.put(readBase, (readCountsMinus.get(readBase) + 1));
						}
						else
						{
							readCountsMinus.put(readBase, 1);
						}
					}

					// Count the strands seen //

					if(strandsSeen.containsKey(readBase))
					{
						String alreadySeen = strandsSeen.get(readBase);
						if(!(alreadySeen.length() >= 2 || alreadySeen.equals(strand)))
						{
							strandsSeen.put(readBase, (strandsSeen.get(readBase) + strand));
						}
					}
					else
					{
						strandsSeen.put(readBase, strand);
					}

					if(qualitySum.containsKey(readBase))
					{
						qualitySum.put(readBase, (qualitySum.get(readBase) + baseQuality));
						mapQualitySum.put(readBase, (mapQualitySum.get(readBase) + mapQuality));
					}
					else
					{
						qualitySum.put(readBase, baseQuality);
						mapQualitySum.put(readBase, mapQuality);
					}
				}
				else
				{
					// Base did not meet quality //
//					System.err.println("Low quality base: " + readBase + " " + baseQuality);
				}

				j++;
				readStart = false;
			}
			// INDEL Processing //
			else if(readBase.equals("+") || readBase.equals("-"))
			{
				String indelType = "";

				if(readBase.equals("+"))
				{
    				indelType = "INS";
				}
				else
				{
    				indelType = "DEL";
				}

				// If the previous base was a reference, count this read as reference but with indel //

				if(prevBase.equals(".") || prevBase.equals(","))
				{
					if(prevBaseQuality >= minAvgQual)
						reads1indel++;
				}

				// Get deletion size and bases //
				int indel_size = 0;
				int max_parse = 1;
				String indelBases = "";
				try {
					String stringWithSize = arrBases[i + 1] + arrBases[i + 2] + arrBases[i + 3];
					stringWithSize = stringWithSize.replaceAll("[^0-9]", "");
					indel_size = Integer.parseInt(stringWithSize);
					max_parse = indel_size + Integer.toString(indel_size).length();

    				for(int bases_parsed = 0; bases_parsed < max_parse; bases_parsed++)
    				{
    					String thisBase = arrBases[i + 1 + bases_parsed];
    					try {
    						Integer.parseInt(thisBase);	// Try to parse an integer from this string, which would be part of indel size
    					}
    					catch (Exception e)
    					{
    						// If no integer, count it.
    						if(thisBase.equals(".") || thisBase.equals(","))
    							bases_parsed = max_parse;
    						else if (thisBase.toUpperCase().equals("A") || thisBase.toUpperCase().equals("C") || thisBase.toUpperCase().equals("G") || thisBase.toUpperCase().equals("T") || thisBase.toUpperCase().equals("N"))
    							indelBases += thisBase;
    					}
    					//indelBases += arrBases[i + 3 + bases_parsed];
    				}
    				// Adjust i to beyond this indel //
    				i = i + max_parse;
				}
				catch (Exception e)
				{
					indel_size = Integer.parseInt(arrBases[i + 1]);
    				for(int bases_parsed = 0; bases_parsed < indel_size; bases_parsed++)
    				{
    					indelBases += arrBases[i + 2 + bases_parsed];
    				}
    				// Adjust i to beyond this indel //
    				i = i + 1 + indel_size;
				}

				// Determine strand //
				if(indelBases.equals(indelBases.toUpperCase()))
				{
					strand = "+";
				}
				else
				{
					strand = "-";
				}

				// Correct case of alleles //
				indelBases = indelBases.toUpperCase();

				// Build an indel key //

				String indelKey = indelType + "-" + indel_size + "-" + indelBases;

				// Count the read //
				if(readCounts.containsKey(indelKey))
				{
					readCounts.put(indelKey, (readCounts.get(indelKey) + 1));
				}
				else
				{
					readCounts.put(indelKey, 1);
				}

				// Count strand-based read count //
				if(strand.equals("+"))
				{
					// Plus strand //
					if(readCountsPlus.containsKey(indelKey))
					{
						readCountsPlus.put(indelKey, (readCountsPlus.get(indelKey) + 1));
					}
					else
					{
						readCountsPlus.put(indelKey, 1);
					}
				}
				else
				{
					// Minus Strand //
					if(readCountsMinus.containsKey(indelKey))
					{
						readCountsMinus.put(indelKey, (readCountsMinus.get(indelKey) + 1));
					}
					else
					{
						readCountsMinus.put(indelKey, 1);
					}
				}

				// Count the strands seen //

				if(strandsSeen.containsKey(indelKey))
				{
					String alreadySeen = strandsSeen.get(indelKey);
					if(!(alreadySeen.length() >= 2 || alreadySeen.equals(strand)))
					{
						strandsSeen.put(indelKey, (strandsSeen.get(indelKey) + strand));
					}
				}
				else
				{
					strandsSeen.put(indelKey, strand);
				}

				if(j < arrQualities.length)
				{
					baseQuality = arrQualities[j] - 33;
					j++;
				}
				if(j < mapQualities.length)
					mapQuality = mapQualities[j] - 33;

				if(qualitySum.containsKey(indelKey))
				{
					qualitySum.put(indelKey, (qualitySum.get(indelKey) + baseQuality));
					mapQualitySum.put(indelKey, (mapQualitySum.get(indelKey) + mapQuality));
				}
				else
				{
					qualitySum.put(indelKey, baseQuality);
					mapQualitySum.put(indelKey, mapQuality);
				}

				readStart = false;
			}
			else if(readBase.toUpperCase().equals("N"))
			{
				// Ignore the base, but keep moving forward for qualities //
				j++;
			}
			else if(readBase.equals("^"))
			{
				// Read start - skip the next base, which is mapping quality //
				i++;
				readStart = true;
			}
			else if(readBase.equals("$"))
			{
				// End of read //
//				i++;
				readStart = false;

			}
			else
			{
				if(readBase.equals(".") || readBase.equals(","))
				{
					// This is the reference base that precedes an indel. Don't advance quality //
				}
				else
				{
					// Ignore characters like * which indicates a pad //
					j++;
				}

			}
		}

		// Declare results hash //
		HashMap<String, String> results = new HashMap<String, String>();

		// Get ref base read counts //

		int strands1 = 0;
		if(strandsSeen.containsKey("ref"))
			strands1 = strandsSeen.get("ref").length();

		// Get average quality //

		int avgQual1 = 0;
		if(reads1 > 0)
			avgQual1 = qualitySum.get("ref") / reads1;

		// Get average map quality //

		int avgMapQual1 = 0;
		if(reads1 > 0)
			avgMapQual1 = mapQualitySum.get("ref") / reads1;

		// Get strand-specific read counts //
		int reads1plus = 0;
		int reads1minus = 0;
		if(readCountsPlus.containsKey("ref"))
			reads1plus = readCountsPlus.get("ref");
		if(readCountsMinus.containsKey("ref"))
			reads1minus = readCountsMinus.get("ref");

		// Append ref info to read counts //
		if(reads1 < 0)
			reads1 = 0;
		results.put(refBase, reads1 + "\t" + strands1 + "\t" + avgQual1 + "\t" + avgMapQual1 + "\t" + reads1plus + "\t" + reads1minus + "\t" + reads1indel);

		// GO through all possible variant keys //

		String[] variantKeys = (String[]) readCounts.keySet().toArray(new String[0]);
		Arrays.sort(variantKeys);
		for(String key : variantKeys)
		{
			int reads2 =  readCounts.get(key);

			// Get strand-specific read counts //
			int reads2plus = 0;
			int reads2minus = 0;
			if(readCountsPlus.containsKey(key))
				reads2plus = readCountsPlus.get(key);
			if(readCountsMinus.containsKey(key))
				reads2minus = readCountsMinus.get(key);

			// Count number of variant-supporting strands //

			int strands2 = 0;
			if(strandsSeen.containsKey(key))
				strands2 = strandsSeen.get(key).length();

			// Get average quality //

			int avg_qual2 = qualitySum.get(key) / reads2;

			// Get average mapping quality //

			int avg_map_qual2 = mapQualitySum.get(key) / reads2;

			if(reads2 > 0)
			{
//				System.err.println("Saving " + key + ": " + reads2 + "\t" + strands2 + "\t" + avg_qual2 + "\t" + avg_map_qual2 + "\t" + reads2plus + "\t" + reads2minus);
				results.put(key, reads2 + "\t" + strands2 + "\t" + avg_qual2 + "\t" + avg_map_qual2 + "\t" + reads2plus + "\t" + reads2minus);
			}
		}

		return(results);
	}



	/**
	 * Counts the depth of read bases meeting a minimum quality
	 *
	 * @param	refBase		Reference base at this position
	 * @param	readBases	String of read bases from pileup
	 * @param	readQuals	String of read base qualities from pileup
	 * @param	minAvgQual	Integer of minimum required base quality to count a base.
	 * @return	results		HashMap<String, String> of results for each allele
	 */
	static int qualityDepth(String readQuals, int minAvgQual)
	{
		int baseQuality = 0;
		int qualityDepth = 0;

		char[] arrQualities = readQuals.toCharArray();

		// Set quality position offset //
		int j = 0;

		// Go through each base //

		for(j = 0; j < arrQualities.length; j++)
		{
				baseQuality = arrQualities[j] - 33;
				if(baseQuality >= minAvgQual)
				{
					qualityDepth++;
				}
		}

		return(qualityDepth);
	}

	/**
	 * Makes the base call (SNP, indel, or consensus) based on read counts
	 *
	 * @param	refBase		Reference base at this position
	 * @param	readCounts	HashMap of read counts for each base observed
	 * @param	callType	Type of call to make (SNP, indel, or consensus)
	 * @param	minReads2	Minimum number of supporting reads to call a variant
	 * @param	minVarFreq	Minimum observed variant frequency to call a variant
	 * @param	minAvgQual	Integer of minimum required base quality to count a base.
	 * @param	pValueThreshold	Significance threshold below which variants will be called
	 * @return	call		The base call made at this position
	 */
	static String callPosition(String refBase, HashMap<String, String> readCounts, String callType, int minReads2, double minVarFreq, int minAvgQual, double pValueThreshold, double minFreqForHom)
	{
		String callResult = "";
		DecimalFormat df = new DecimalFormat("###.##");

		int reads1 = 0;
		int reads2 = 0;
		int readsWithIndels = 0;
		int strands1 = 0;
		int strands2 = 0;
		int avgQual1 = 0;
		int avgQual2 = 0;
		int avgMap1 = 0;		// Average mapping quality of reference-supporting reads
		int avgMap2 = 0;		// Average mapping quality of variant-supporting reads
		int reads1indel = 0;	// Reference-supporting reads that contain indel at next base
		int reads1plus = 0;		// Reference-supporting reads on plus strand
		int reads1minus = 0;	// Reference-supporting reads on minus strand
		int reads2plus = 0;		// Variant-supporting reads on plus strand
		int reads2minus = 0;	// Variant-supporting reads on minus strand
		double pValue = 1;
		double varFreq = 0.00;
		String varAllele = "";

		try
		{
			if(readCounts.containsKey(refBase))
			{
				try
				{
					String[] refBaseContents = readCounts.get(refBase).split("\t");
					reads1 = Integer.parseInt(refBaseContents[0]);
					strands1 = Integer.parseInt(refBaseContents[1]);
					avgQual1 = Integer.parseInt(refBaseContents[2]);
					avgMap1 = Integer.parseInt(refBaseContents[3]);
					reads1plus = Integer.parseInt(refBaseContents[4]);
					reads1minus = Integer.parseInt(refBaseContents[5]);

					if(refBaseContents.length > 6)
						reads1indel = Integer.parseInt(refBaseContents[6]);
				}
				catch(Exception e)
				{
					System.err.println("Error parsing refBase readcounts from " + readCounts.get(refBase));
				}
			}

			String[] alleleKeys = (String[]) readCounts.keySet().toArray(new String[0]);

			Arrays.sort(alleleKeys);

			// Get the total number of reads at this position //
			int totalReadCounts = 0;
			for(String allele : alleleKeys)
			{
				String[] alleleContents = readCounts.get(allele).split("\t");
				try {
					int thisReads = Integer.parseInt(alleleContents[0]);
					totalReadCounts += thisReads;
				}
				catch(Exception e)
				{
				}
			}

			for(String allele : alleleKeys)
			{
				String[] alleleContents = readCounts.get(allele).split("\t");

				if(allele.equals(refBase))
				{
					// Skip the reference base; we got that already //
				}
				else
				{
					// Reset variables //

					int thisReads1 = reads1;
					int thisReads2 = 0;
					int thisStrands2 = 0;
					int thisAvgQual2 = 0;
					int thisAvgMap2 = 0;
					int thisReads2plus = 0;
					int thisReads2minus = 0;

					// Parse the information //

					try {
						thisReads2 = Integer.parseInt(alleleContents[0]);
						thisStrands2 = Integer.parseInt(alleleContents[1]);
						thisAvgQual2 = Integer.parseInt(alleleContents[2]);
						thisAvgMap2 = Integer.parseInt(alleleContents[3]);
						thisReads2plus = Integer.parseInt(alleleContents[4]);
						thisReads2minus = Integer.parseInt(alleleContents[5]);
						// If this is an indel, make note of it //

						if(allele.contains("INS") || allele.contains("DEL"))
						{
							readsWithIndels += thisReads2;
						}
					}
					catch (Exception e)
					{

					}


					if(!callType.equals("CNS") || thisReads2 > reads2)
					{
						//double thisVarFreq = (double) thisReads2 / (double) (reads1 + thisReads2);
						double thisVarFreq = (double) thisReads2 / (double) totalReadCounts;
						double thisPvalue = 1;
						// For indels, adjust the read1 count //
						if(allele.contains("INS") || allele.contains("DEL"))
						{
							//System.err.println(allele + " gets " + thisReads2 + " " + thisVarFreq);
							// Adjust the reads1 counts which include reads supporting indels //
//							thisReads1 = reads1 - reads1indel;
//							if(thisReads1 < 0)
//								thisReads1 = 0;

//							thisVarFreq = (double) thisReads2 / (double) (thisReads1 + thisReads2);

							// Correct for indel-containing reads, but ensure we don't overcorrect //
							int thisTotalReadCounts = totalReadCounts - reads1indel;
							if(thisTotalReadCounts < thisReads2)
								thisTotalReadCounts = thisReads2;

							// Compute new variant allele frequency //
							thisVarFreq = (double) thisReads2 / (double) totalReadCounts;
						}

						// Calculate the p-value //
						if(pValueThreshold == 0.99)
						{
							thisPvalue = 0.98;
						}
						else
						{
							thisPvalue = getSignificance(reads1, thisReads2);
						}


						// Save the most frequent variant allele, even if we won't use it //
						if(thisReads2 > reads2 && thisAvgQual2 >= minAvgQual)
						{
//							System.err.println(allele + " passed with " + thisReads2);
							if(allele.contains("INS") || allele.contains("DEL"))
							{
								varAllele = getShortIndel(allele);
							}
							else
							{
								varAllele = allele;
							}

							reads2 = thisReads2;
							strands2 = thisStrands2;
							avgQual2 = thisAvgQual2;
							avgMap2 = thisAvgMap2;
							reads2plus = thisReads2plus;
							reads2minus = thisReads2minus;
							varFreq = thisVarFreq * 100;
							pValue = thisPvalue;
						}
						else
						{
							//System.err.println(allele + " failed with " + thisReads2 + " " + thisAvgQual2);
						}

						// Call the variant if it meets calling criteria //

						if(thisReads2 >= minReads2 && thisAvgQual2 >= minAvgQual && thisVarFreq >= minVarFreq)
						{
							thisReads1 = reads1;
							thisVarFreq = thisVarFreq * 100;

							// Determine type of variant //
							String thisVarType = "SNP";
							if(allele.contains("INS") || allele.contains("DEL"))
							{
								thisVarType = "INDEL";
								thisReads1 = reads1;
								if(thisReads1 < 0)
									thisReads1 = 0;
								// Change allele to short indel version //
								allele = getShortIndel(allele);
							}

							if(thisPvalue <= pValueThreshold)
							{
								// Call the variant if we're variant calling //
								if(callType.equals("SNP") || callType.equals("INDEL"))
								{
									reads2 = thisReads2;
									strands2 = thisStrands2;
									avgQual2 = thisAvgQual2;
									avgMap2 = thisAvgMap2;
									reads2plus = thisReads2plus;
									reads2minus = thisReads2minus;
									pValue = thisPvalue;

									// Convert to consensus-like genotype //

									String genotype = "";
									if(thisVarFreq >= (minFreqForHom * 100))
									{
										genotype = allele + allele;
										if(thisVarType.equals("INDEL"))
											genotype = allele + "/" + allele;
									}
									else
									{
										genotype = refBase + allele;
										if(thisVarType.equals("INDEL"))
											genotype = "*/" + allele;
									}

									// Only report the desired variant type //

									if(thisVarType.equals(callType))
									{
										// Report the variant regardless //
										if(callResult.length() > 0)
											callResult += "\n";

										if(thisReads1 < 0)
											thisReads1 = 0;

										if(reads2 < 0)
											reads2 = 0;

										//callResult += allele + "\t" + reads1 + "\t" + reads2 + "\t" + df.format(thisVarFreq) + "%\t" + strands1 + "\t" + strands2 + "\t" + avgQual1 + "\t" + avgQual2 + "\t" + pValue;
										callResult += genotypeToCode(genotype) + "\t" + thisReads1 + "\t" + reads2 + "\t" + df.format(thisVarFreq) + "%\t" + strands1 + "\t" + strands2 + "\t" + avgQual1 + "\t" + avgQual2 + "\t" + pValue;
										callResult += "\t" + avgMap1 + "\t" + avgMap2;
										callResult += "\t" + reads1plus + "\t" + reads1minus + "\t" + reads2plus + "\t" + reads2minus + "\t" + varAllele;
									}

								}
								else if(callType.equals("CNS") && thisReads2 >= reads2)
								{
									reads2 = thisReads2;
									strands2 = thisStrands2;
									avgQual2 = thisAvgQual2;
									avgMap2 = thisAvgMap2;
									reads2plus = thisReads2plus;
									reads2minus = thisReads2minus;
									pValue = thisPvalue;

									String genotype = "";
									if(thisVarFreq >= (minFreqForHom * 100))
									{
										genotype = allele + allele;
										if(thisVarType.equals("INDEL"))
											genotype = allele + "/" + allele;
									}
									else
									{
										genotype = refBase + allele;
										if(thisVarType.equals("INDEL"))
											genotype = "*/" + allele;
									}

									callResult = genotypeToCode(genotype) + "\t" + thisReads1 + "\t" + reads2 + "\t" + df.format(thisVarFreq) + "%\t" + strands1 + "\t" + strands2 + "\t" + avgQual1 + "\t" + avgQual2 + "\t" + pValue;
									callResult += "\t" + avgMap1 + "\t" + avgMap2;
									callResult += "\t" + reads1plus + "\t" + reads1minus + "\t" + reads2plus + "\t" + reads2minus + "\t" + varAllele;
								}

							}
							else
							{
								// Somehow p-value not less than threshold //

							}
						}
						else
						{
							// Did not meet reads2, variant allele frequency, base quality, or p-value thresholds //
						}

					}



				}
			}
		}
		catch(Exception e)
		{
	    	System.err.println("Read Counts Exception: " + e.getLocalizedMessage());
	    	e.printStackTrace(System.err);
		}

		// If we must have a call result for CNS calling, decide on reference or NO call //
		if(callResult.length() == 0 && callType.equals("CNS"))
		{

			if(reads1 > 0 && reads1 > minReads2)
			{
				// Call reference because enough reads supporting ref base were observed //
				callResult = refBase + "\t" + reads1 + "\t" + reads2 + "\t" + df.format(varFreq) + "%\t" + strands1 + "\t" + strands2 + "\t" + avgQual1 + "\t" + avgQual2 + "\t" + pValue;
				callResult += "\t" + avgMap1 + "\t" + avgMap2;
				callResult += "\t" + reads1plus + "\t" + reads1minus + "\t" + reads2plus + "\t" + reads2minus + "\t" + varAllele;
			}
			else
			{
				callResult = "N" + "\t" + reads1 + "\t" + reads2 + "\t" + df.format(varFreq) + "%\t" + strands1 + "\t" + strands2 + "\t" + avgQual1 + "\t" + avgQual2 + "\t" + pValue;
				callResult += "\t" + avgMap1 + "\t" + avgMap2;
				callResult += "\t" + reads1plus + "\t" + reads1minus + "\t" + reads2plus + "\t" + reads2minus + "\t" + varAllele;
			}

		}

		return(callResult);
	}


	/**
	 * Evaluates strand bias in variant allele read counts according to provided p-value threshold
	 *
	 * @param	reads1plus	Number of reference-supporting reads on + strand
	 * @param	reads1minus	Number of reference-supporting reads on + strand
	 * @param	reads2plus	Number of reference-supporting reads on + strand
	 * @param	reads2minus	Number of reference-supporting reads on + strand
	 * @param	strandPvalueThreshold	P-value threshold below which variant fails strand filter
	 * @return	call		A string with strand filter status, counts, p-value
	 */
	static String strandFilter(int reads1plus, int reads1minus, int reads2plus, int reads2minus, double strandPvalueThreshold)
	{
		DecimalFormat pvalueFormat = new DecimalFormat("0.####E0");
		String strandFilterStatus = "Pass:" + reads1plus + ":" + reads1minus + ":" + reads2plus + ":" + reads2minus + ":" + 1;
		double refStrandPlus = 0.50;
		double varStrandPlus = 0.50;
		double strandPvalue = 1.00;

		// Calculate strandedness for variant allele //

		if((reads2plus + reads2minus) > 0)
			varStrandPlus = (double) reads2plus / (double) (reads2plus + reads2minus);

		// To save time, only calculate p-value if var strandedness is biased //

		if(varStrandPlus < 0.10 || varStrandPlus > 0.90)
		{
			// Calculate strandedness for reference allele if we have 2+ reads //

			if((reads1plus + reads1minus) > 1)
			{
				refStrandPlus = (double) reads1plus / (double) (reads1plus + reads1minus);
				strandPvalue = VarScan.getSignificance(reads1plus, reads1minus, reads2plus, reads2minus);
			}
			// Otherwise, only homozygous-variant reads seen, so compare to a 50/50 distribution //
			else
			{
				// Compare to expected 50/50 distribution //
				int testReads1plus = (int) (reads2plus + reads2minus) / 2;
				int testReads1minus = (reads2plus + reads2minus) - testReads1plus;
				strandPvalue = VarScan.getSignificance(testReads1plus, testReads1minus, reads2plus, reads2minus);
			}

			strandFilterStatus = "Pass:" + varStrandPlus + ":" + reads1plus + ":" + reads1minus + ":" + reads2plus + ":" + reads2minus + ":" + pvalueFormat.format(strandPvalue);

			// If ref allele had good strandedness, and var allele did not, this may be a failure //
			if(refStrandPlus >= 0.10 && refStrandPlus <= 0.90 && !(varStrandPlus >= 0.10 && varStrandPlus <= 0.90))
			{
				if(strandPvalue < strandPvalueThreshold)
				{
					strandFilterStatus = "Fail:" + reads1plus + ":" + reads1minus + ":" + reads2plus + ":" + reads2minus + ":" + pvalueFormat.format(strandPvalue);
				}
			}
		}

		return(strandFilterStatus);
	}


	/**
	 * Calculates significance of read counts versus baseline error
	 *
	 * @param	obsReads1	Reads supporting allele 1
	 * @param	obsReads2	Reads supporting allele 2
	 * @return	p-value	P-value from Fisher's Exact Test
	 */
	public static double getSignificance(int obsReads1, int obsReads2)
	{
		double pValue = 1;
		double baseline_error = 0.01;

		int coverage = obsReads1 + obsReads2;

		int expReads2 = (int) (coverage * baseline_error);
		int expReads1 = coverage - expReads2;

		pValue = getSignificance(expReads1, expReads2, obsReads1, obsReads2);
		return(pValue);
	}


	/**
	 * Calculates significance of read counts between two samples
	 *
	 * @param	expReads1	Reads supporting allele 1 (expected)
	 * @param	expReads2	Reads supporting allele 2 (expected)
	 * @param	obsReads1	Reads supporting allele 1 (observed)
	 * @param	obsReads2	Reads supporting allele 2 (observed)
	 * @return	p-value 	P-value from Fisher's Exact Test
	 */
	public static double getSignificance(int expReads1, int expReads2, int obsReads1, int obsReads2)
	{
		double pValue = 1;

		if(expReads1 < 0)
			expReads1 = 0;

		if(expReads2 < 0)
			expReads2 = 0;

		if(obsReads1 < 0)
			obsReads1 = 0;

		if(obsReads2 < 0)
			obsReads2 = 0;

		// Set up fisher's exact test //

		FishersExact fisher = new FishersExact(expReads1 + expReads2 + obsReads1 + obsReads2 + 100);

		// Calculate a p-value //

		pValue = fisher.getRightTailedP(expReads1, expReads2, obsReads1, obsReads2);
		int fisher_max = 1000;
		int num_tries = 0;

		while(Double.isNaN(pValue) && num_tries < 10)
		{
			fisher = new FishersExact(expReads1 + expReads2 + obsReads1 + obsReads2 + fisher_max);
			//pValue = fisher.getTwoTailedP(expReads1, expReads2, obsReads1, obsReads2);
			pValue = fisher.getRightTailedP(expReads1, expReads2, obsReads1, obsReads2);
			fisher_max = fisher_max + 1000;
			num_tries++;
		}

		if(num_tries >= 10)
			System.err.println("Warning: unable to calculate p-value failure: " + expReads1 + "," + expReads2 + "," + obsReads1 + "," + obsReads2);

		// If p-value is 1, do left-sided test //

		if(pValue >= 0.999)
		{
			pValue = fisher.getLeftTailedP(expReads1, expReads2, obsReads1, obsReads2);

			while(Double.isNaN(pValue))
			{
				fisher = new FishersExact(expReads1 + expReads2 + obsReads1 + obsReads2 + fisher_max);
				//pValue = fisher.getTwoTailedP(expReads1, expReads2, obsReads1, obsReads2);
				pValue = fisher.getLeftTailedP(expReads1, expReads2, obsReads1, obsReads2);
				fisher_max = fisher_max + 1000;
			}
		}

		return(pValue);
	}


	/**
	 * Converts from two-allele genotype to IUPAC code
	 *
	 * @param	genotype
	 * @return	code    	P-value from Fisher's Exact Test
	 */
	static String genotypeToCode(String gt)
	 {
		 if(gt.equals("AA"))
			 return("A");

		 if(gt.equals("CC"))
			 return("C");

		 if(gt.equals("GG"))
			 return("G");

		 if(gt.equals("TT"))
			 return("T");

		 if(gt.equals("AC") || gt.equals("CA"))
			 return("M");

		 if(gt.equals("AG") || gt.equals("GA"))
			 return("R");

		 if(gt.equals("AT") || gt.equals("TA"))
			 return("W");

		 if(gt.equals("CG") || gt.equals("GC"))
			 return("S");

		 if(gt.equals("CT") || gt.equals("TC"))
			 return("Y");

		 if(gt.equals("GT") || gt.equals("TG"))
			 return("K");

		 if(gt.substring(0, 1).equals("N"))
			 return("N");

		 return(gt);
	 }


	/**
	 * Gets variant allele from ref Base and consensus call
	 *
	 * @param	genotype
	 * @return	code    	P-value from Fisher's Exact Test
	 */
	static String getVarAllele(String refBase, String consCode)
	 {
		String varAllele = consCode;

		if(consCode.contains("/"))
		{
			String[] tempArray = consCode.split("/");
			if(tempArray.length > 1)
				varAllele = tempArray[1];
		}
		else if(consCode.equals("M") || consCode.equals("R") || consCode.equals("W") || consCode.equals("S") || consCode.equals("Y") || consCode.equals("K"))
		{
			if(consCode.equals("M"))
			{
				if(refBase.equals("A"))
					varAllele = "C";
				else
					varAllele = "A";
			}

			else if(consCode.equals("R"))
			{
				if(refBase.equals("A"))
					varAllele = "G";
				else
					varAllele = "A";
			}

			else if(consCode.equals("W"))
			{
				if(refBase.equals("A"))
					varAllele = "T";
				else
					varAllele = "A";
			}

			else if(consCode.equals("S"))
			{
				if(refBase.equals("G"))
					varAllele = "C";
				else
					varAllele = "G";
			}

			else if(consCode.equals("Y"))
			{
				if(refBase.equals("C"))
					varAllele = "T";
				else
					varAllele = "C";
			}

			else if(consCode.equals("K"))
			{
				if(refBase.equals("G"))
					varAllele = "T";
				else
					varAllele = "G";
			}
		}

		return(varAllele);
	 }

	/**
	 * Converts from long to short indel allele
	 *
	 * @param	genotype
	 * @return	code    	P-value from Fisher's Exact Test
	 */
	static String getShortIndel(String gt)
	{
		 if(gt.contains("INS") || gt.contains("DEL"))
		 {
			 try
			 {
				 String[] gtContents = gt.split("-");
				 String indel_type = gtContents[0];
				 String indel_size = gtContents[1];
				 String indel_bases = gtContents[2];

				 if(indel_type.contains("INS"))
				 {
					 return("+" + indel_bases);
				 }
				 else
				 {
					 return("-" + indel_bases);
				 }
			 }
			 catch(Exception e)
			 {
				 System.err.println("Warning: error generating consensus from " + gt);
			 }

		 }

		 return("N");
	 }


	/**
	 * Returns true if a variant is heterozygous
	 *
	 * @param	genotype   		An ambiguity code or a a1/a2 genotype
	 * @return	boolean			True if heterozygous
	 */
	 static boolean isHeterozygous(String genotype)
	 {
		 if(genotype.contains("/"))
		 {
			 String[] alleles = genotype.split("/");
			 if(!alleles[0].equals(alleles[1]))
				 return(true);
		 }
		 else
		 {
			 if(genotype.equals("M") || genotype.equals("R") || genotype.equals("W") || genotype.equals("S") || genotype.equals("Y") || genotype.equals("K"))
				 return(true);
		 }

		 return(false);
	 }

	/**
	 * Returns true if a variant is homozygous
	 *
	 * @param	genotype   		An ambiguity code or a a1/a2 genotype
	 * @return	boolean			True if homozygous
	 */
	 static boolean isHomozygous(String genotype)
	 {
		 if(genotype.contains("/"))
		 {
			 String[] alleles = genotype.split("/");
			 if(alleles[0].equals(alleles[1]))
				 return(true);
		 }
		 else
		 {
			 if(genotype.equals("A") || genotype.equals("C") || genotype.equals("G") || genotype.equals("T"))
				 return(true);
		 }

		 return(false);
	 }


	/**
	 * Converts IUPAC codes to two-letter genotypes
	 *
	 * @param	genotype   		An ambiguity code or a a1/a2 genotype
	 * @return	boolean			True if homozygous
	 */
	 static String codeToGenotype(String code)
	 {
			if(code.equals("A") || code.equals("C") || code.equals("G") || code.equals("T"))
				return(code + "/" + code);
			else if(code.contains("/"))
			{
				return(code);
			}
			else if(code.equals("M"))
				return("A/C");
			else if(code.equals("R"))
				return("A/G");
			else if(code.equals("W"))
				return("A/T");
			else if(code.equals("S"))
				return("C/G");
			else if(code.equals("Y"))
				return("C/T");
			else if(code.equals("K"))
				return("G/T");
			else
			return("N/N");
	 }



}



