
import sys
from sys import argv
import subprocess
import os
import argparse 
from datetime import datetime 

def pipeline(sampleID, reference, ranking_method, consensus_bul, *optionalfilters): #  must put the optional filters in string format
	
	sample = sampleID 
	firstsample = sample 
	starting_index = False 
	ending_index = False
	percentage_thresh = False 
	pval_thresh = False 
	read2_thresh = False 
	ranking_method = str(ranking_method)

	for argument in optionalfilters: # arguments needed for the ranking script  
		if "-s" in argument: 
			starting_index = argument.strip("-s")
		if "-e" in argument:  
                        ending_index = argument.strip("-e")
	 	if "-perc" in argument:  
                        percentage_thresh = argument.strip("-perc")
		if "-pval" in argument:  
                        pval_thresh = argument.strip("-pval")
		if "-r" in argument:  
                        read2_thresh = argument.strip("-r")        
	
	# these initial two inputs are being read from the input directory. All the next iteration inputs will come from this directory itself
	if consensus_bul is False: 
		h5file = "../input/" + sampleID + ".fastq" 
		filereference= "../input/" + reference
	if consensus_bul is True:
		h5file =  sampleID + ".fastq" # this is for the heading consensus files that are created in the output directory
		filereference = reference  
	
	
	print "THIS IS SAMPLE ID " + sampleID 
	#filereference  = reference 	
	samfile =  str(sample) + ".sam"
    	unalignedfile = str(sample) + ".unaligned.fasta"
    	bamfile =   str(sample) + ".sorted.bam"
    	sortfile =   str(sample) + ".sorted"
   	m5file =   str(sample) + ".m5"
    	pileupfile =   str(sample) + ".pileup"
    	indelfile =   str(sample)+ ".indel"
        readcountsfile =  str(sample) + ".readcounts"
    	snpfile = str(readcountsfile) +  ".snp"
	consensusfile =  str(sample).strip(".realign") + ".consensus.fasta"
        rankedfile = str(sample) + ".rankedsnp"
		
		# first realignment to reference 
	
	if consensus_bul is False: 
		command = "bwa index " + filereference
	command = "blasr " + h5file + " " + filereference + " -minMatch 19 -nproc 16 -bestn 1 -unaligned "+ unalignedfile+ " -m 5 -out " +m5file
        print command
        subprocess.call(command, shell = True)

        command = "/opt/Sparc/compiled/Sparc b "+filereference+" m "+m5file+" k 2 g 2 c 2 t 0.1 o "+sample
        print command
        print >> sys.stderr
        subprocess.call(command, shell =True)	
	
	sample = sampleID + ".realign"
	#h5file =  sampleID  + ".fastq"
	samfile =  str(sample) + ".sam"
    	unalignedfile = str(sample) + ".unaligned.fasta"
    	bamfile =   str(sample) + ".sorted.bam"
    	sortfile =   str(sample) + ".sorted"
   	m5file =   str(sample) + ".m5"
   	pileupfile =   str(sample) + ".pileup"
   	indelfile =   str(sample)+ ".indel"
   	readcountsfile =  str(sample) + ".readcounts"
   	snpfile = str(readcountsfile) +  ".snp"
	consensusfile =  str(sampleID) + ".consensus.fasta"
   	rankedfile = str(sample) + ".rankedsnp"
	
	#  command = "bwa index " + consensusfile
	
	print "THIS IS H5 file" + h5file  + " " + consensusfile
	
	command = "blasr " + h5file + " " + consensusfile + " -nproc 10 -bestn 1 -unaligned "+ unalignedfile+ " -sam -out " +samfile
	print command
    	subprocess.call(command, shell = True)
    
    	#command = "bwa mem -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0 -t 30 " + consensusfile + " " + h5file + " |samtools view -Sb - |samtools sort - -o " + bamfile
	#creates sam|bam|sortedbam
	command = "samtools view -Sb " + samfile+ " |samtools sort - -o " + bamfile 
	print command	
	print >> sys.stderr 
	subprocess.call(command, shell =True)
	print >> sys.stderr, "\n bam_file is done!\n" + str(bamfile) +  "\n"
		
	command = "rm " + samfile
        print command
        print >> sys.stderr
        subprocess.call(command, shell =True)

	command = "samtools index" + " " + bamfile # keep 
	print command	
	print >> sys.stderr
	subprocess.call(command, shell = True)
	
	command = "samtools mpileup -BQ0 -d100000 -f"+" "+ consensusfile  + " " + bamfile +" > " + pileupfile
	print command 
	print >> sys.stderr
	subprocess.call(command, shell = True)
	print >> sys.stderr,"\n pileup_file is done!\n" + str(pileupfile) + "\n"
	
	command = "java -jar ../script/VarScan.v2.4.2.jar readcounts" + " "+ pileupfile + " " + " --output-file" + " " + readcountsfile + " --min-base-qual 5"
	print command	
	print >> sys.stderr
	subprocess.call(command, shell = True)
	print >> sys.stderr, "\n readcounts_file done!\n" + str(readcountsfile)  +   "\n"
	
	command = "java -cp ../script/ ReadCount_consensus " + readcountsfile # running as background process. This creates the consensusfile. We actually declared it before, but now it has actually been filled  
	print command
	subprocess.call(command, shell = True)
	print >> sys.stderr,"\n consensus file is done!\n" + consensusfile + "\n"
	print "java clas DONE!"        
		

	command  = "java -cp ../script/ ReadCount_snp" + " " + readcountsfile
	print command	
	print >> sys.stderr 
	subprocess.call(command, shell = True)
	print >> sys.stderr, "\n snp_file is done!\n"  + str(snpfile) + "\n"


	command = " python rankSNVs.py" + " " + str(snpfile)+ " " + rankedfile  + " " + ranking_method + " " + consensusfile + " "  #" " + "-s" + " " + str(args.starting_index) + " " + "-e" + " " + str(args.ending_index) + " " + "-perc" + " " + str(args.percentage) + " " + "-pval" + " " + str(args.pval_thresh) + " " + "-r "+ " " + str(args.read2_thresh) 
	#output file is where the ranked snips are going to be, along with the bam file, this will feed into the cluster program 
	
	if starting_index :
        	command = command + "-s" + " " + starting_index  + " "
	if ending_index :
        	command = command + "-e" + " " + ending_index + " "
	if percentage_thresh:
        	command = command + "-perc" + " " + percentage_thresh+ " "
	if pval_thresh:
        	command = command + "-pval" + " " + pval_thresh + " "
	if read2_thresh:
        	command = command + "-r" + " " + read2_thresh + " "
   	print command
	
	subprocess.call(command, shell = True)


	print "FINSIHED"
	endtime = datetime.now() 
#	print endtime - starttime 

	returnlist = [readcountsfile, snpfile, consensusfile, bamfile, rankedfile, sample, firstsample] 
	return returnlist  




	#return second_alignment  # note the consensusfile that is being passed on to the next subcluster is only the initial reference for the clusters that came before, we are still going to have to realign to the consensus of that specific subcluster 
# returns  ranked, bam, consensus, snp => To be used in clustering program and next iteration 
