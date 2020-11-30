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


start = time.time()
print "Post processing consensus sequences"
outDir="/mnt/LHRI/Bioinformatics/Projects/Quasispecies/lhri_quasi-seq/testing/19mix_Miseq/hivalign/"
dataDir=outDir
#outDir="/mnt/LHRI/Bioinformatics/Projects/Quasispecies/lhri_quasi-seq/testing/"
#dataDir="/mnt/LHRI/Bioinformatics/Projects/Quasispecies/manuscript/docker/lhri_quasi-seq/quasiSeqOut/5b79ae5b-f902-41fe-95cd-60df8b54cdda/"
#outDir=dataDir
finalConsensusFile=dataDir+"Vpr.fasta"
#finalConsensusFile=dataDir+"final_consensus.fasta"
comparisonFile=outDir+"post_processed_consensus_comparison.txt"
finalIdFile=outDir+"post_processed_consensus_ids.txt"
threshold=99.5
print outDir
print finalConsensusFile
command = "blasr "+finalConsensusFile+" "+finalConsensusFile+" -m 1 | sort -u -k2 | awk '{ split( $1, a, \"/\" );if (a[1]!=$2 && $6>99.9) print a[1],$2,$6}' > " + comparisonFile
print "Post-processing blasr command: "+command
subprocess.call(command, shell = True,executable="/bin/bash")
print >> sys.stderr, "\n Post-processing blasr is done!\n" + str(outDir) + "final_consensus_comparison.txt\n"
comparisonSeqs={}
with open(comparisonFile, 'r') as f:
	for line in f:
		comparison = line.split()
		if float(comparison[2])>=threshold:
			currSet=frozenset([comparison[0],comparison[1]])
			if comparison[0] in comparisonSeqs:
				comparisonSeqs[comparison[0]]=frozenset(currSet | comparisonSeqs[comparison[0]])
			else:
				comparisonSeqs[comparison[0]]=currSet
			if comparison[1] in comparisonSeqs:
				comparisonSeqs[comparison[1]]=frozenset(currSet | comparisonSeqs[comparison[1]])
			else:
				comparisonSeqs[comparison[1]]=currSet
encountered=[]
finalFqRecords=[]
finalIds=[]
with open(finalConsensusFile, 'r') as f:
	for record in SeqIO.parse(f,'fasta'):
		if not record.id in encountered:
		   currFqRecords=[]
		   if record.id in comparisonSeqs:
			   finalIds.append('~'.join(map(str, comparisonSeqs[record.id])))
			   for id in comparisonSeqs[record.id]:
				   encountered.append(id)
				   '''with open(dataDir+id, 'r') as fq:
					   for fqRecord in SeqIO.parse(fq,'fastq'):
						   currFqRecords.append(fqRecord)'''
		   else:
			   encountered.append(record.id)
			   finalIds.append(record.id)
			   '''with open(dataDir+record.id, 'r') as fq:
				   for fqRecord in SeqIO.parse(fq,'fastq'):
					   currFqRecords.append(fqRecord)'''

		   #finalFqRecords.append(currFqRecords)
with open(finalIdFile, 'w') as output:
	for consensus_id in finalIds:
		currFqRecords=[]
		output.write(consensus_id+"\n");
		ids=consensus_id.split('~')
		for id in ids:
			with open(dataDir+record.id, 'r') as fq:
				for fqRecord in SeqIO.parse(fq,'fastq'):
					currFqRecords.append(fqRecord)


end = time.time()
print "writeFinalConsensusStats took " + str((end - start))+ "seconds"
