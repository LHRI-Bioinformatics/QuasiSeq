# source /home/software/pitchfork/deployment/setup-env.sh
from __future__ import division 
import gevent.monkey
gevent.monkey.patch_all()
import gevent
import gevent.pool

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed 
from sparc import pipeline
import argparse
from multiprocessing import Pool, freeze_support, Manager
import subprocess
import os
import time 
from collections import defaultdict 
import multiprocessing.pool 

numtime = 1
starttime = time.time()
sampleQueue = defaultdict(str) 
#pool = Pool() 

#current_directory = os.getcwd()
parser = argparse.ArgumentParser()
parser.add_argument("filename" , help = "input file needed" )
parser.add_argument("reference", help = "initial reference file needed")
#parser.add_argument("outputfile",help = "output file needed" )
parser.add_argument("ranking_method", choices = ["pval", "percentage"], default = "percentage", help = "first preference ranking method")

# these next arguments are all optional filters . There are default values if any of these are not specified by the user
parser.add_argument("-s","--starting_index",type = int, default = 0, help = "optional starting index")
parser.add_argument("-e","--ending_index", type = int, default = None, help = "optional ending index")
parser.add_argument("-perc","--percentage_thresh", type = float, default = 1.0, help = "optional lowest threshold for percentage")
parser.add_argument("-pval","--pval_thresh", type = float, default = .001, help = "optional highest threshold for pvalue")
parser.add_argument("-r","--read2_thresh", type = int, default = 10, help = "optional lowest threshold for the second read")
parser.add_argument("-snp", "--snp_thresh", type = int, default = 10, help = "The min requirement of ranked SNP's needed for a subcluster iteration")
args = parser.parse_args()

filename = str(args.filename)
reference = str(args.reference)
#output = str(args.outputfile)
ranking = str(args.ranking_method)
snpthresh = args.snp_thresh

if args.starting_index is not 0:
        startingindex = str(args.starting_index) + "-s"
else :
        startingindex = str(args.starting_index)
if args.ending_index is not None :
        endingindex = str(args.ending_index) + "-e"
else :
        endingindex = str(args.ending_index)
if args.percentage_thresh is not 1.0:
        percentagethresh = str(args.percentage_thresh) + "-perc"
else:
        percentagethresh = str(args.percentage_thresh)
if args.pval_thresh is not .001 :
        pvalthresh =str(args.pval_thresh) + "-pval"
else:
        pvalthresh = str(args.pval_thresh)
if args.read2_thresh is not 10 :
        read2thresh =str(args.read2_thresh) + "-r"
else :
        read2thresh = str(args.read2_thresh)
	

def mediary(tuple, pool): 
	print "mediary started"

	new_Input = tuple[0] 
	reference = tuple[1]
	pool = pool 
	print "This is new_Input" + new_Input  
	print "This is reference" + reference 
	print "calling this mediary function "

	fullrecursionpipeline(new_Input, reference, ranking, startingindex, endingindex, percentagethresh, pvalthresh, read2thresh, snpthresh, True, pool)
	
	return 



def fullrecursionpipeline(*args):
	
	 global starttime 
         global numtime
         filename = args[0]
         reference = args[1]
         ranking = args[2]
         startingindex = args[3]
         endingindex = args[4]
         percentagethresh = args[5]
         pvalthresh = args[6]
         read2thresh = args[7]
	 snpthresh = args[8]

	 conbul = args[9]
	 pool = args[10]
	 global sampleQueue  
	  
         cluster_parameters = pipeline(filename, reference, ranking,conbul, startingindex, endingindex, percentagethresh, pvalthresh, read2thresh)  #the initial run will have False for the consenus_bul because we are using the original 

         conbul = True

     	 rankedfile = str(cluster_parameters[4].strip())

    	 
 	 with open(rankedfile  ,"r") as ranked:
		for line in ranked: 
			ranked_amount = int(line.split(None, 1)[0])
			print "THIS IS NUMBER OF RANKED " + str(ranked_amount)
  			break 
	 
	 if ranked_amount  >= snpthresh:
		adjustedQueue = dict(sampleQueue)  
	 	for key in adjustedQueue: 
	 		if key == str(filename):
	 			del sampleQueue[filename]
	 			
		command = "java -classpath .:/usr/local/MATLAB/MATLAB_Runtime/v91/toolbox/javabuilder/jar/javabuilder.jar:/opt/quasi-seq/multi_driver/script/SigClust.jar getClust /opt/quasi-seq/multi_driver/data/" + cluster_parameters[3] + " /opt/quasi-seq/multi_driver/data/" + cluster_parameters[4] + " " + cluster_parameters[6]  
		print command
		subprocess.call(command, shell =True)
		cluster_directory = str(cluster_parameters[6]) + "_clusters"
		print ("THIS IS CLUSTER DIRECTORY: " + cluster_directory)
		
 
	 	with open(cluster_directory,"r") as clusteroutput: 
				
			for line in clusteroutput: 
				clusterheading = line.strip() 
				sampleQueue[clusterheading] = False 
		with open(cluster_directory, "r") as clusteroutput: 
			for line in clusteroutput:
				clusterlist = [] 
                                clusterheading = line.strip()
                                command = "java -cp ../script/ Fastq_filter_by_readID " + filename+ ".fastq " + clusterheading # the cluster program will create clusterheading files in the output folder so thats what you are using here. Output/clusterheading is the actual file that is located in the output folder. remember the clusterheading b$
                                print command
                                print "analyzing the " + str(numtime) + "SUBCLUSTERs"
                                subprocess.call(command, shell = True)
                                newInput =   clusterheading
                                print "THIS IS THE NEXT INPUT" + newInput
                                newReference = cluster_parameters[2] # the consensus file of the previous cluster
                                numtime +=1
                                newTuple = (newInput, newReference )
                                clusterlist.append(newTuple)
                                print "************************Recursive Call*************************"
				pool.spawn(mediary, newTuple, pool)
				
		print str(time.time()-starttime)
          
	 else:			
		sampleQueue[cluster_parameters[6]] = True # the initial sample is done, has no more subclusters in it
		print "THIS IS SAMPLE QUEUE"
		print cluster_parameters[6]
		print sampleQueue 
		
		command = "cat " + str(cluster_parameters[2]) + " >> " + str(consensus_string)
		subprocess.call(command, shell = True)
		
		command = "echo " + str(cluster_parameters[6]) + " >> " + str(headings_string) 
		
		consensus_dict[cluster_parameters[6]]  =cluster_parameters[2]  
 
 		
		subprocess.call(command, shell = True)



		
		bul_counter = 0 
		for value in sampleQueue.values() :
			#print "checking"
			if value == True:
				bul_counter+=1
			#	print "BUL COUNT" 
			#	print bul_counter 
			#	print len(sampleQueue)
				if bul_counter ==  len(sampleQueue):
					print "END OF PIPELINE " + finalconsensus_string + " " + finalheadings_string 
					print "Pipeline took " + str(time.time()-starttime)+ " seconds"
					total_lines = 0 
					for key in sampleQueue: 
						total_lines += sum(1 for line in open(key))
					
					
					for key in sampleQueue:
						print key
						headings_output.write(str(key))
						headings_output.write("\n")

						num_lines = sum(1 for line in open(key))	
					
						for element in consensus_dict : 
							if key == element: 
								#consensus_output.write(str(element) + "\n")
								with open(str(consensus_dict[element]), "r") as consensus:
									#next(consensus)
									for line in consensus:
										for c in line: 
											if c is ">":
												perc_int = (num_lines/total_lines) * 100
												perc_int = round(perc_int, 2)
												#print "this is num_lines " + str(num_lines)
												#print "this is total lines " + str(total_lines) 
												num_string = " Freq: " + str(perc_int) + "% "
												#print "this is nums_string" + num_string
												line = line.replace(".consensus", num_string)
												#line += str(num_lines)
												#lines = line.strip("\n") 
										consensus_output.write(line)
										
						
						
		
					print sampleQueue
                                        #command = "rm " + finalheadings_string + " " + finalconsensus_string
                                        #subprocess.call(command, shell = True) 
					command = "python matcher.py " + headings_string + " " + consensus_string
                                        subprocess.call(command, shell = True)


if __name__ == "__main__":
		
	
	var = subprocess.Popen("nproc", stdout = subprocess.PIPE) 
	nproc = var.stdout.read() # number of cores on machine
	print "total processors: "
	print nproc
        #if nproc < 16: 
        	#nproc = 16
        print "MAIN FUNCTION STARTS..."        
        consensus_string = filename + "_Quasi_strains.fasta"
        finalconsensus_string = filename + "_final" + "_Quasi_strains.fasta"
        
        clearing_string = open((filename + "_Quasi_strains.fasta"), "w")
        clearing_string.close() 

	clearing_consensusstring = open(finalconsensus_string, "w")
	clearing_consensusstring.close()  
        
        headings_string = filename + "_leaves"
        finalheadings_string = filename + "_final" + "_leaves"

	clearing_headingsstring = open(headings_string, "w")
	clearing_headingsstring.close() 
        
        consensus_output = open(finalconsensus_string, "a+")
	headings_output = open(finalheadings_string, "w+")
	
	consensus_dict = defaultdict(str)
  	pool = gevent.pool.Pool(int(nproc))  
        fullrecursionpipeline(filename, reference, ranking, startingindex, endingindex, percentagethresh, pvalthresh, read2thresh, snpthresh, False, pool)
	pool.join() 
