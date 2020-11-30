import logging
import sys, getopt, time
from sys import argv
import subprocess
import os
import argparse
from datetime import datetime
import errno
import csv
from operator import itemgetter

import pysam
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from shutil import copyfile

import numpy
import re

gene="Vif"

start = time.time()
print "Post processing consensus sequences"
outDir="/mnt/LHRI/Bioinformatics/Projects/Quasispecies/lhri_quasi-seq/testing/19mix_Miseq/hivalign/"
dataDir=outDir
#outDir="/mnt/LHRI/Bioinformatics/Projects/Quasispecies/lhri_quasi-seq/testing/"
#dataDir="/mnt/LHRI/Bioinformatics/Projects/Quasispecies/manuscript/docker/lhri_quasi-seq/quasiSeqOut/5b79ae5b-f902-41fe-95cd-60df8b54cdda/"
#outDir=dataDir

finalConsensusFile=dataDir+str(gene)+".fasta"

#finalConsensusFile=dataDir+"final_consensus.fasta"
comparisonFile=outDir+str(gene)+"_comparison.txt"
matrixFile=outDir+str(gene)+"_matrix.txt"
threshold=99.5
print outDir
print finalConsensusFile

command = "blasr "+finalConsensusFile+" "+finalConsensusFile+" -m 1 | sort -u -k2 | awk '{ split( $1, a, \"/\" );print a[1],$2,$6}' > " + comparisonFile

print "Post-processing blasr command: "+command
subprocess.call(command, shell = True,executable="/bin/bash")
print >> sys.stderr, "\n Post-processing blasr is done!\n" + str(comparisonFile)+"\n"
matrix={}
encountered=set()
with open(comparisonFile, 'r') as f:
	for line in f:
		comparison = line.split()
		if comparison[0] not in encountered:
			encountered.add(comparison[0])
			matrix[comparison[0]]={}
		if comparison[1] not in encountered:
			encountered.add(comparison[1])
			matrix[comparison[1]]={}
			
		matrix[comparison[0]][comparison[1]]=comparison[2]
		matrix[comparison[1]][comparison[0]]=comparison[2]
		
with open(matrixFile, 'w') as f:
	#header
	for id1 in encountered:
		f.write("\t"+str(id1))
	f.write("\n")
	
	#data
	for id1 in encountered:
		f.write(str(id1))
		for id2 in encountered:
			if id1 in matrix:
				#print "Found "+str(id1)
				if id2 in matrix[id1]:
					f.write("\t"+str(matrix[id1][id2]))
					#print "Found "+str(id2) + " in "+str(id1)
		   			
		   		else:
		   			f.write("\tNA")
		   			#print "Could not find "+str(id2) + " in "+str(id1)
		   			
		   	#else:
		   		#print "Could not find "+str(id1)
		f.write("\n")   		



end = time.time()
print "writeFinalConsensusStats took " + str((end - start))+ "seconds"
