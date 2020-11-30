


from simplified_module import pipeline
import argparse
from multiprocessing import Pool, freeze_support, Manager
import subprocess
import os
import time 
#from concurrent import futures 
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
parser.add_argument("outputfile",help = "output file needed" )
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
output = str(args.outputfile)
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


		
class NoDaemonProcess(multiprocessing.Process):
	  # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

class MyPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess
	

def mediary(tuple): 
	print tuple 
	filename = tuple[0] 
	reference = tuple[1]
	print "This is filename" + filename 
	print "This is reference" + reference 
	print "calling this mediary function "
	passed_dictionary = tuple[2]
	fullrecursionpipeline(filename, reference, ranking, startingindex, endingindex, percentagethresh, pvalthresh, read2thresh, snpthresh, passed_dictionary, True)
	
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
	 passed_dictionary = args[9]
	 conbul = args[10]
	# myPool = args[11]
	 clusterlist = [] 
	 global sampleQueue  
	  
         cluster_parameters = pipeline(filename, reference, ranking,conbul, startingindex, endingindex, percentagethresh, pvalthresh, read2thresh)  #the initial run will have False for the consenus_bul because we are using the original 

         # you need to run the cluster program here with cluster_parameters, only the first and second elements, as the parameters # the cluster program is going to have to return the path of the file that it creates which we will store as cluster_directory 
         conbul = True
         # proceed only if the cluster directory file is not empty 
     	 rankedfile = str(cluster_parameters[4].strip())
     	 print "rankedfile:", rankedfile
     	 
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
	 			
		command = "java -classpath .:/usr/local/MATLAB/MATLAB_Runtime/v91/toolbox/javabuilder/jar/javabuilder.jar:/mnt/LHRI/Bioinformatics/Projects/Quasispecies/manuscript/scripts/SigClust/for_testing/SigClust.jar getClust /mnt/LHRI/Bioinformatics/Projects/Quasispecies/manuscript/rishub_backup/xiaoli/" + cluster_parameters[3] + " /mnt/LHRI/Bioinformatics/Projects/Quasispecies/manuscript/rishub_backup/xiaoli/" + cluster_parameters[4] + " " + cluster_parameters[6]  
		print command
		subprocess.call(command, shell =True)
		cluster_directory = str(cluster_parameters[6]) + "_clusters"
		print ("THIS IS CLUSTER DIRECTORY: " + cluster_directory)
		
		#if os.path.exists(cluster_directory):'
		#print "TEST"     
	 	with open(cluster_directory,"r") as clusteroutput: #where cluster output is the output of the cluster function
				
			for line in clusteroutput: 
				print "TEST"	
				clusterheading = line.strip() 
				sampleQueue[clusterheading] = False 
		with open(cluster_directory, "r") as clusteroutput: 
			for line in clusteroutput:

                                print "TEST 3"
                                clusterheading = line.strip()
                                print clusterheading
                                command = "java -cp ../script/ Fastq_filter_by_readID ../input/10clones_FL_redirected.fastq " + clusterheading # the cluster program will create clusterheading files in the output folder so thats what you are using here. Output/clusterheading is the actual file that is located in the output folder. remember the clusterheading b$
                                print command
                                print "analyzing the" + str(numtime) + "SUBCLUSTER"
                                subprocess.call(command, shell = True)
                                newInput =   clusterheading
                                print "THIS IS THE NEXT INPUT" + newInput
                                newReference = cluster_parameters[2] # the consensus file of the previous cluster
                                numtime +=1
                                print "THIS IS THE " + str(numtime) + "th subcluster"
                                #print "This cluster took" + str(time.time()-starttime) 
                                #passed_dictionary[newInput] = False
                                newTuple = (newInput, newReference, passed_dictionary, )
                                clusterlist.append(newTuple)
                                print "************************Recursive Call*************************"
                                        #testTuple = (newInput, newReference, ranking, startingindex, endingindex, percentagethresh, pvalthresh, read2thresh, snpthresh, passed_dictionary, True, myPool)
                                        #pool.map(fullrecursionpipeline, testTuple)
                                fullrecursionpipeline(newInput, newReference, ranking, startingindex, endingindex, percentagethresh, pvalthresh, read2thresh, snpthresh, passed_dictionary, True)

	 		print "TEST 2"
                	numclusters = len(clusteroutput.readlines())
			print numclusters
			#if numclusters > 0 :
				#print "TEST 3" 
			numtime = 0
			
			for lin4 in clusteroutput:
				print "TEST 4"			
			#for line in clusteroutput:				
			#	print "TEST 3"				
			#	clusterheading = line.strip()
			#	print clusterheading
                        #	command = "java -cp ../script/ Fastq_filter_by_readID ../input/10clones_FL_redirected.fastq " + clusterheading # the cluster program will create clusterheading files in the output folder so thats what you are using here. Output/clusterheading is the actual file that is located in the output folder. remember the clusterheading by itself is just the name of the file you need the path to get to the actual file
                        #	print command 
			#	print "analyzing the" + str(numtime) + "SUBCLUSTER"
			#	subprocess.call(command, shell = True)
			#	newInput =   clusterheading
		        #	print "THIS IS THE NEXT INPUT" + newInput
		        #	newReference = cluster_parameters[2] # the consensus file of the previous cluster
			#	numtime +=1 
			#	print "THIS IS THE " + str(numtime) + "th subcluster"
			#	#print "This cluster took" + str(time.time()-starttime) 
			#	passed_dictionary[newInput] = False 
			#	newTuple = (newInput, newReference, passed_dictionary, )
			#	clusterlist.append(newTuple)
			#	print "************************Recursive Call*************************"
					#testTuple = (newInput, newReference, ranking, startingindex, endingindex, percentagethresh, pvalthresh, read2thresh, snpthresh, passed_dictionary, True, myPool)
					#pool.map(fullrecursionpipeline, testTuple)
		        #	fullrecursionpipeline(newInput, newReference, ranking, startingindex, endingindex, percentagethresh, pvalthresh, read2thresh, snpthresh, passed_dictionary, True)
		        	
				#myPool.map(mediary, clusterlist)
			
			

		print str(time.time()-starttime)
                       	# print "THIS IS" + cluster_parameters[2]
                       	# print "THIS IS " + str(numtime) + "running"
                       	# p.map(fullrecursionpipeline(newInput, newReference, newOutput, ranking, startingindex, endingindex, percentagethresh, pvalthresh, read2thresh), [])
                       	# print "RECURSING THROUGH SUB SUB CLUSTER"
                       	#p.map(fullrecursionpipeline, newInput, newReference, newOutput)
              
          	#p.map(fullrecursionpipeline, newInput, newReference, newOutput, ranking, conbul, startingindex, endingindex, percentagethresh, pvalthresh, read2thresh)

   		#else:
                #	print "FAILED sample queue"
		#	passed_dictionary[cluster_parameters[6]] = True # the initial sample is done, has no more subclusters in it
                #	print "THIS IS"
                #	print cluster_parameters[6]
                #	print passed_dictionary.values()

                #	command = "cat " + str(cluster_parameters[2]) + " >> " + str(consensus_string)
                #	subprocess.call(command, shell = True)
#
 #               	headings_output.write(str(cluster_parameters[6]) + "\n")
  #                                      #headings_output.close()
#
  #              	bul_counter = 0
 # #              	for value in passed_dictionary.values() :
    #                    	if value == True:
     #                   	        bul_counter+=1
      #                          	if bul_counter ==  len(passed_dictionary):
       #                                 	print "END OF PIPELINE " + consensus_string + " " + headings_string
	
	 else:			
		sampleQueue[cluster_parameters[6]] = True # the initial sample is done, has no more subclusters in it
		print "THIS IS SAMPLE QUEUE"
		print cluster_parameters[6]
		print sampleQueue 
		
		command = "cat " + str(cluster_parameters[2]) + " >> " + str(consensus_string)
		subprocess.call(command, shell = True)

		command = "cat" + str(cluster_parameters[6]) + " >> " + str(headings_string) 
 
		subprocess.call(command, shell = True)


		
		#headings_output.write(str(cluster_parameters[6]) + "\n")
					#headings_output.close()
		
		bul_counter = 0 
		for value in sampleQueue.values() :
			print "checking"
			if value == True:
				bul_counter+=1
				print "BUL COUNT" 
				print bul_counter 
				print len(sampleQueue)
				if bul_counter ==  len(sampleQueue):
					print "END OF PIPELINE " + consensus_string + " " + headings_string 
					
					
					
					for key in sampleQueue:
						print key
						headings_output.write(str(key))
						headings_output.write("\n")
						
 					#headings_output.close() 
						
		
					print sampleQueue 

if __name__ == "__main__":
        print "THIS IS MAIN FUNCTION"
        manager = Manager() 
        freeze_support()
        passed_dictionary = manager.dict() 
        consensus_string = filename + "_Quasi_strains.fasta"
        
        clearing_string = open((filename + "_Quasi_strains.fasta"), "w")
        clearing_string.close()  
        
        headings_string = filename + "_leavesnew"
        #consensus_output = open(consensus_string, "a+")
        headings_output = open(headings_string, "w+")
  #    	pool = Pool() 
        fullrecursionpipeline(filename, reference, ranking, startingindex, endingindex, percentagethresh, pvalthresh, read2thresh, snpthresh, passed_dictionary, False)
		
