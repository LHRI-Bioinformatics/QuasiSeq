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

import multiprocessing as mp

def getFullLengthReads(outDir,sampleName,m5File, start, end, buffer):
	saved_args = locals()
    	print("getFullLengthReads parameters:", saved_args)
	m5 = open(m5File,"r")
	fullLengthReads=open(outDir+sampleName+"_fullLength","w")
	nonFullLengthReads=open(outDir+sampleName+"_nonFullLength","w")
	print "start-end:"+str(start)+"-"+str(end)
	startBuffer=start+buffer
	endBuffer=end-buffer
	# print "startBuffer:endBuffer:"+str(startBuffer)+"-"+str(endBuffer)
	for line in m5:
		data = line.split()
		readName=data[0]

		# if len(data)<9:
		# 	print(readName)

		if int(data[7])<=(startBuffer) and int(data[8])>=endBuffer:
			fullLengthReads.write(readName+'\t'+data[9]+'\n')
		else:
			nonFullLengthReads.write(readName+'\t'+data[9]+'\n')
	fullLengthReads.close()
	nonFullLengthReads.close()

def getFinalSparcClusterConsensus(outDir,clusterName,polishedConsensusFile):
	referenceFile=outDir+clusterName+'.consensus.fasta'
	fastq=outDir+clusterName+'.fastq'
	m5File=outDir+clusterName+".m5"
	#finalRoundConsensusFile=outDir+clusterName+'.finalRoundSparc.consensus.fasta'
	alignFastqBlasrM5(fastq,referenceFile,m5File)# Align subreads to the initial reference and index resulting bam
   	getConsensusSparc(m5File,referenceFile,polishedConsensusFile)# Get consensus using Sparc (first round)

def truncateFasta(fastaIn, fastaOut, start, end):
	truncated_sequences = [] # Setup an empty list
	for record in SeqIO.parse(fastaIn, "fasta"):
	    truncated_sequences.append(SeqRecord(record.seq[start:end],record.id, '', ''))

	SeqIO.write(truncated_sequences, fastaOut, "fasta")

def writeFinalConsensusStats(outDir):
	start = time.time()
	print "$$$$$$$$$$$$$$$$$$###################!!!!!!!!!!Writing final consensus stats"
	readStatsFile = outDir+'final_consensus_stats.txt'
	statsExist=os.path.isfile(readStatsFile)
	data=list()
	total = 0

	if statsExist:
		with open(readStatsFile, 'rb') as f:
			reader = csv.reader(f)
			for line in reader:
				data.append(line)
				total += int(line[0])

		data.sort(key=lambda x: int(x[0]),reverse=True)  # 1 being the column number
		totalRecord = [total, 'TOTAL']
		data.append(totalRecord)

		with open(readStatsFile, 'wb') as f:
    			csv.writer(f,lineterminator='\n').writerows(data)

def writeFinalClusterConsensusStats(outDir,clusterName):
	start = time.time()
	print "!!!!!!!!!!Writing final consensus stats for "+clusterName
	readStatsFile = outDir+'final_consensus_stats.txt'
	with open(readStatsFile, 'a+') as finalConsensusStatsFile:
			finalConsensusStatsFile.write(str(len((open(outDir+clusterName,"r")).readlines()))+','+clusterName+'\n')

def writeFinalClusterConsensus(outDir,clusterName):
	start = time.time()
	print "!!!!!!!!!!Writing final consensus for "+clusterName

	#polishedConsensusFile=outDir+clusterName+'.firstRoundSparc.consensus.fasta'
	#getFinalSparcClusterConsensus(outDir,clusterName,polishedConsensusFile)

	polishedConsensusFile=outDir+clusterName+'.finalRoundSparc.consensus.fasta'
	getFinalSparcClusterConsensus(outDir,clusterName,polishedConsensusFile)

	clusterConsensus = SeqIO.read(polishedConsensusFile, "fasta")
	finalConsensusFile=open(outDir+'final_consensus.fasta',"a")
	finalConsensusFile.write('>'+clusterName+'\n')
	finalConsensusFile.write(str(clusterConsensus.seq)+'\n')
	writeFinalClusterConsensusStats(outDir,clusterName)

def alignFastqAndIndexSam_v1_3(fastq,reference,bamFile, *optionalfilters):
	start = time.time()
	outFile=os.path.splitext(bamFile)[0]
	samFile=outFile+".sam"

	command = "blasr " + fastq + " " + reference + " -nproc 30 -bestn 1 -sam -out " +samFile

	if writeUnaligned:
		unalignedOut = outFile + ".unaligned.fa"
		command = command + " -unaligned " + unalignedOut

	print command
	print >> sys.stderr
	subprocess.call(command, shell = True,executable="/bin/bash")
	print >> sys.stderr, "\n sam_file is done!\n" + str(samFile) +  "\n"

	#command = "/opt/bbmap/reformat.sh sam=1.3 in="+ samFile + " out=stdout.sam -Xmx20g | samtools view -Sb - | samtools sort - -o " + bamFile
	command = "samtools view -Sb "+ samFile + " | samtools sort - -o " + bamFile
	print command
	print >> sys.stderr
	subprocess.call(command, shell = True,executable="/bin/bash")
	print >> sys.stderr, "\n bam_file is done!\n" + str(bamFile) +  "\n"

	#command = "rm " + samFile
	#print command
	#print >> sys.stderr
	#subprocess.call(command, shell =True)

	command = "samtools index" + " " + bamFile # keep
	print command
	print >> sys.stderr
	subprocess.call(command, shell = True,executable="/bin/bash")

	end = time.time()
	print "alignFastqAndIndexPacBioSam_v1_3 took " + str((end - start))+ "seconds"

def alignFastqBlasrM5(fastq,reference,m5File, *optionalfilters):
	start = time.time()
	command = "blasr " + fastq + " " + reference + " -minMatch 19 -nproc 30 -bestn 1 -m 5 -out " +m5File

	if writeUnaligned:
		unalignedOut = m5File + ".ref_unaligned.fa"
		command = command + " -unaligned " + unalignedOut

	print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"+command
	print >> sys.stderr
	subprocess.call(command, shell = True,executable="/bin/bash")
	print >> sys.stderr, "\n m5File is done!\n" + str(m5File) +  "\n"

	end = time.time()
	print "alignFastqBlasrM5 took " + str((end - start))+ "seconds"

def normalizeStartEnd(consensusFile,referenceFile, startIndex, endIndex, *optionalfilters):
	####################Bad logic; do not use!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	start = time.time()

	m1File=os.path.splitext(consensusFile)[0]+".m1"
	command = "blasr " + consensusFile + " " + referenceFile + " -minMatch 19 -nproc 30 -bestn 1 -m 1 -out " +m1File

	if writeUnaligned:
		unalignedOut = m1File + ".ref_unaligned.fa"
		command = command + " -unaligned " + unalignedOut

	print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"+command
	print >> sys.stderr
	subprocess.call(command, shell = True,executable="/bin/bash")
	print >> sys.stderr, "\n m1File is done!\n" + str(m1File) +  "\n"

	with open(m1File,"r") as m1:
		data = m1.readline().split()
	end = time.time()
	print "normalizeStartEnd took " + str((end - start))+ "seconds"
	return [(int(startIndex)+(int(data[9])-int(data[6]))),(int(endIndex)+(int(data[10])-int(data[7])))]

def indexReference(reference, *optionalfilters):

	start = time.time()
	command = "bwa index " + reference
	print command
	#subprocess.call(command, shell = True)
	end = time.time()

	start = time.time()
	command = "samtools faidx " + reference
	print command
	subprocess.call(command, shell = True)
	end = time.time()

	print "indexReference took " + str((end - start))+ "seconds"


def getClusterFastq(fastq,currCluster,currClusterFastq, *optionalfilters):

	start = time.time()
	#saved_args = locals()
    	#print("getClusterFastq args:", saved_args)
    	#clusterReadNames=set(open(currCluster, 'r').read().splitlines())

	clusterReadNames = set()
	regexp = re.compile(r'/[0-9]+_[0-9]+')
	with open(currCluster, 'r') as f:
		for line in f:
			clean_line=line.rstrip('\n')
			if regexp.search(clean_line):
				clusterReadNames.add(clean_line[:clean_line.rindex('/')])
			else:
				clusterReadNames.add(clean_line)


	#print(clusterReadNames)
	print(str(len(clusterReadNames))+" reads in "+currCluster)
	with open(fastq) as allReads, open(currClusterFastq, 'w') as clusterReads:
	    records = SeqIO.parse(allReads, 'fastq')
	    for record in records:
		#print record.id
		if str(record.id) in clusterReadNames:
			#print("Found "+record.id)
			SeqIO.write(record, clusterReads, 'fastq')

	end = time.time()
	print "getClusterFastq took " + str((end - start))+ "seconds"

def getOrientedFastq(fastq,readNames,readNamesFastq, *optionalfilters):

	start = time.time()
	#saved_args = locals()
    	#print("getClusterFastq args:", saved_args)
    	clusterReadNames = {}
    	regexp = re.compile(r'/[0-9]+_[0-9]+')
	with open(readNames, 'r') as f:
		for line in f:
			(key, val) = line.split()
			if regexp.search(key):
				clusterReadNames[key[:key.rindex('/')]] = val
			else:
				clusterReadNames[key] = val
    	#clusterReadNames=set(open(readNames, 'r').read().splitlines())

	#print(clusterReadNames)
	print(str(len(clusterReadNames))+" reads in "+readNames)
	with open(fastq) as allReads, open(readNamesFastq, 'w') as clusterReads:
	    records = SeqIO.parse(allReads, 'fastq')
	    for record in records:
		#print record.id
		if str(record.id) in clusterReadNames:
			#print("Found "+record.id)
			if clusterReadNames[str(record.id)]=='+':
				SeqIO.write(record, clusterReads, 'fastq')
			else:
				#SeqIO.write(record.reverse_complement(id=record.id, description=str(record.description)+' reverse_complement'), clusterReads, 'fastq')
				SeqIO.write(record.reverse_complement(id=record.id, description=record.description), clusterReads, 'fastq')

	end = time.time()
	print "getOrientedFastq took " + str((end - start))+ "seconds"


def getConsensusLength(fasta, *optionalfilters):

	start = time.time()
	consensus = SeqIO.read(fasta, "fasta")
	end = time.time()
	print "getConsensusLength took " + str((end - start))+ "seconds"
	return len(consensus.seq)

def getSNPs(readCountsFile, estimatedErrorRate, readCountSNPPvalThresh, readCountSNPVarFreqThresh, readCountSNPCoverageThresh, *optionalfilters):

	start = time.time()
	snpFile=readCountsFile+".snp"  #os.path.splitext(readCountsFile)[0]+".snp"
	command = "java -cp ./ ReadCount_snp " + readCountsFile + " " +snpFile+ " " +str(estimatedErrorRate)+ " " +str(readCountSNPPvalThresh)+ " " +str(readCountSNPCoverageThresh)+ " " +str(readCountSNPVarFreqThresh)
	# check the snpfile , if em
	print(command)
	print >> sys.stderr
	subprocess.call(command, shell = True)

	end = time.time()
	print "getSNPs took " + str((end - start))+ "seconds"

def getReadCounts(reference,bamFile, estimatedErrorRate, readCountSNPPvalThresh, readCountSNPVarFreqThresh, readCountSNPCoverageThresh, *optionalfilters):

	start = time.time()
	outFile=os.path.splitext(bamFile)[0]
	pileupFile=outFile+ ".pileup"
	readCountsFile=outFile+ ".readcounts"



	command = "samtools mpileup -BQ0 -d100000 -f"+" "+ reference  + " " + bamFile +" > " + pileupFile
	print command
	print >> sys.stderr
	subprocess.call(command, shell = True)
	print >> sys.stderr,"\n pileup_file is done!\n" + str(pileupFile) + "\n"

	command = "java -jar VarScan.v2.4.2.jar readcounts" + " "+ pileupFile + " " + " --output-file" + " " + readCountsFile + " --min-base-qual 5"
	print command
	print >> sys.stderr
	subprocess.call(command, shell = True)
	print >> sys.stderr, "\n readcounts_file done!\n" + str(readCountsFile)  +   "\n"
	end = time.time()
	print "getReadCounts took " + str((end - start))+ "seconds"

	getSNPs(readCountsFile, estimatedErrorRate, readCountSNPPvalThresh, readCountSNPVarFreqThresh, readCountSNPCoverageThresh)

def getPileup(reference,bamFile, *optionalfilters):

	start = time.time()
	outFile=os.path.splitext(bamFile)[0]
	pileupFile=outFile+ ".pileup"
	readCountsFile=outFile+ ".readcounts"

	command = "samtools mpileup -BQ0 -d100000 -f"+" "+ reference  + " " + bamFile +" > " + pileupFile
	print command
	print >> sys.stderr
	subprocess.call(command, shell = True)
	print >> sys.stderr,"\n pileup_file is done!\n" + str(pileupFile) + "\n"
	end = time.time()
	print "getPileup took " + str((end - start))+ "seconds"

def getConsensusSparc(m5File, referenceFile, consensusFile, *optionalfilters):

	start = time.time()
	#command = "/opt/Sparc/compiled/Sparc b "+referenceFile+" m "+m5File+" k 2 g 2 c 2 t 0.1 o "+consensusFile.replace('.consensus.fasta','')
	command = "/opt/Sparc/compiled/Sparc b "+referenceFile+" m "+m5File+" k 1 g 1 c 2 t 0.2 o "+consensusFile.replace('.consensus.fasta','')
	print command
	subprocess.call(command, shell = True)
	print >> sys.stderr,"\n consensus file is done!\n" + consensusFile + "\n"
	print "Created Consensus!"
	#removeConsensusNs(consensusFile)
	indexReference(consensusFile) # index the consensus file

	end = time.time()
	print "getConsensusSparc took " + str((end - start))+ "seconds"

def consenseAndRealignSparc(fastq,referenceFile,initialBamfile, readCountsFile, consensusFile, realignedBamfile, estimatedErrorRate, readCountSNPPvalThresh, readCountSNPVarFreqThresh, readCountSNPCoverageThresh, *optionalfilters):

	start = time.time()

	print("consenseAndRealignSparc:"+fastq+", "+referenceFile+", "+initialBamfile+", "+readCountsFile+", "+consensusFile+", "+realignedBamfile)

	m5File=str(os.path.splitext(initialBamfile)[0])+".m5"

	#final sparc consensus
   	alignFastqBlasrM5(fastq,referenceFile,m5File)# Align subreads to the first round Sparc consensus and index resulting bam

   	getConsensusSparc(m5File,referenceFile,consensusFile)# Get consensus using Sparc (second round)

	alignFastqAndIndexSam_v1_3(fastq,consensusFile,realignedBamfile)

	getReadCounts(consensusFile,realignedBamfile, estimatedErrorRate, readCountSNPPvalThresh, readCountSNPVarFreqThresh, readCountSNPCoverageThresh) # Get SNPs file for realignment using Xiaoli and Rishub's code
	getPileup(consensusFile,realignedBamfile) # Get pileup from bam

	end = time.time()
	print "consenseAndAlign took " + str((end - start))+ "seconds"


def rankSNVs(snpFile, rankedFile, rankingMethod, consensusFile, startIndex, endIndex, percentThresh, pvalThresh, minorAlleleReadThresh, maxSignatures, *optionalfilters):

	saved_args = locals()
    	print("rankSNVs parameters:", saved_args)

	#snpFile=(os.path.split(snpFile))[1]

	#command = " python rankSNVs.py" + " " + snpFile + " " + rankedFile + " " + rankingMethod + " " + consensusFile +" -s " + str(startIndex) + " " + "-e" + " " + str(endIndex) + " " + "-perc" + " " + str(percentThresh) + " " + "-pval" + " " + str(pvalThresh) + " -r "+ str(minorAlleleReadThresh)
	#print "~~~~~~~~~~~~~~~~~~"+command
	#output file is where the ranked snps are going to be, along with the bam file, this will feed into the cluster program

	command = " python rankSNVs.py" + " " + str(snpFile)+ " " + rankedFile  + " " + rankingMethod + " " + consensusFile + " "


	#output file is where the ranked snips are going to be, along with the bam file, this will feed into the cluster program
	#09112017 startIndex=1
	#09112017 endIndex=getConsensusLength(consensusFile)
	if startIndex :
		command = command + "-s" + " " + str(startIndex)  + " "
	if endIndex :
		command = command + "-e" + " " + str(endIndex) + " "
	if percentThresh:
		command = command + "-perc" + " " + str(percentThresh)+ " "
	if pvalThresh:
		command = command + "-pval" + " " + str(pvalThresh) + " "
	if minorAlleleReadThresh:
		command = command + "-r" + " " + str(minorAlleleReadThresh) + " "

	command = command + "-numSig" + " " + str(maxSignatures) + " "

	print("~~~~~~~~~~~~~~~~~~"+command)
	subprocess.call(command, shell = True)

def rankSNVs_2(pileupFile, rankedFile, rankingMethod, consensusFile, startIndex, endIndex, percentThresh, pvalThresh, minorAlleleReadThresh, *optionalfilters):

	outFile=pileupFile + "2cns.txt"
	command = "java -jar VarScan.v2.4.2.jar pileup2snp " + pileupFile + " --variants --min-avg-qual 5"
	if percentThresh:
		command = command + " --min-var-freq " + str(percentThresh)
	if pvalThresh:
		command = command + " --p-value " + str(pvalThresh)
	if minorAlleleReadThresh:
		command = command + " --min-reads2 " + str(minorAlleleReadThresh)
	command=command+" > "+ outFile

	print ("!!!!######$$$$$$!!!!!" +command)
	print >> sys.stderr
	subprocess.call(command, shell = True)
	print >> sys.stderr, "\n pileup2snp done!\n" + str(outFile)  +   "\n"

	pileup2snp=open(outFile, 'r').readlines() #list(open(outFile, 'r').read().splitlines())
	rankedSNPs=open(rankedFile, 'w')
	rankedSNPs.write(str(len(pileup2snp)-1) + " # of ranked snps with reference "+consensusFile+"\n")
	lastSNP=len(pileup2snp)-2
	if lastSNP>450 : lastSNP=450
	#for line in sorted(pileup2snp[1:], key=lambda line: line.split()[6], reverse=True):
	#for line in sorted(pileup2snp[1:lastSNP], key=lambda line: float(line.split()[11])):
	for line in sorted(pileup2snp[1:], key=lambda line: float(line.split()[11])):
		line = line.strip()
		columns = line.split()
		snpPosition = columns[1]
		if int(snpPosition)<=endIndex:
			rankedSNPs.write(snpPosition+"\n")

def rankSNVs_3(pileupFile, rankedFile, rankingMethod, consensusFile, startIndex, endIndex, percentThresh, pvalThresh, minorAlleleReadThresh, *optionalfilters):

	buffer=30
	outFile=pileupFile + "2snp.txt"
	command = "java -jar VarScan.v2.4.2.jar pileup2snp " + pileupFile + " --variants --min-avg-qual 5"
	if percentThresh:
		command = command + " --min-var-freq " + str(percentThresh)
	if pvalThresh:
		command = command + " --p-value " + str(pvalThresh)
		print("*&^%$#Setting pvalThresh to "+str(pvalThresh)+"*&^%$#")
	#if minorAlleleReadThresh:
	#	command = command + " --min-reads2 " + str(minorAlleleReadThresh)
	command=command+" > "+ outFile

	print ("!!!!######$$$$$$!!!!!" +command)
	print >> sys.stderr
	subprocess.call(command, shell = True)
	print >> sys.stderr, "\n pileup2snp done!\n" + str(outFile)  +   "\n"

	pileup2snp=open(outFile, 'r').readlines() #list(open(outFile, 'r').read().splitlines())

	rankedSNPs=[]

	totalReads=[]
	consensusLength=getConsensusLength(consensusFile)
	for line in pileup2snp[1:]:
		line = line.strip()
		columns = line.split()
		totalReads.append(int(columns[4])+int(columns[5]))
	medianReads=float(numpy.median(sorted(totalReads)))

	#for line in sorted(pileup2snp[1:], key=lambda line: line.split()[6], reverse=True):
	#for line in sorted(pileup2snp[1:], key=lambda line: float(line.split()[11])):

	for line in sorted(pileup2snp[1:], key=lambda line: float((line.split()[6]).replace('%','')), reverse=True):
		line = line.strip()
		columns = line.split()
		snpPosition = columns[1]
		numReads1=int(columns[4])
		numReads2=int(columns[5])
		lineReads=int(columns[4])+int(columns[5])
		#if (lineReads/medianReads)>0.99:
		#	print("("+columns[4]+"+"+columns[5]+")/"+str(medianReads)+"="+str(lineReads/medianReads)+" is greater than 0.99")
		#	if int(snpPosition) not in rankedSNPs and int(snpPosition)<endIndex:
		#		rankedSNPs.append(int(snpPosition))
		if numReads1<int(minorAlleleReadThresh) or numReads2<int(minorAlleleReadThresh) :
			print "******Read1/Read2 logic:"+snpPosition+":"+str(numReads1)+":"+str(numReads2)
		#if int(snpPosition) not in rankedSNPs and int(snpPosition)<=(endIndex-30) and int(snpPosition)>=(startIndex+30) and numReads1>=int(minorAlleleReadThresh) and numReads2>=int(minorAlleleReadThresh):
		if int(snpPosition) not in rankedSNPs and int(snpPosition)<=(consensusLength-buffer) and int(snpPosition)>=buffer:# and numReads1>=int(minorAlleleReadThresh) and numReads2>=int(minorAlleleReadThresh):
			rankedSNPs.append(int(snpPosition))

	rankedSNPsFile=open(rankedFile, 'w')
	rankedSNPsFile.write(str(len(rankedSNPs)) + " # of ranked snps with reference "+consensusFile+"\n")
	for snp in rankedSNPs:
  		rankedSNPsFile.write("%i\n" % snp)

def cluster(outDir,realignedBamfile, rankedFile, sampleName, snpThresh, *optionalfilters):

	with open(rankedFile, 'r') as f:
		first_line = f.readline()
	numSNPs = int(first_line.split(None, 1)[0])
	if not snpThresh: snpThresh=10
	cwd = ""
	print "numSNPs:"+str(numSNPs)+">= snpThresh:"+str(snpThresh)+"->"+str(numSNPs >= snpThresh)
	if numSNPs >= snpThresh:
		print str(numSNPs)+" SNPs: Continuing with clustering"
		command = "java -classpath .:/usr/local/MATLAB/MATLAB_Runtime/v91/toolbox/javabuilder/jar/javabuilder.jar:SigClust.jar getClust " + realignedBamfile + " " + rankedFile + " " + outDir+sampleName
		print command
		subprocess.call(command, shell = True,executable="/bin/bash")
		return True
	else:
		print str(numSNPs)+" SNPs: Creating consensus for "+sampleName+ "(cluster)"
		#command = "cat "+sampleName+".consensus.fasta >> final_consensus.fasta"
		#command = 'cat <(echo ">'+sampleName+'") <(tail -n\+2 '+sampleName+'.consensus.fasta) >> final_consensus.fasta'
		#subprocess.call(command, shell = True,executable="/bin/bash")
		writeFinalClusterConsensus(outDir,sampleName)
		#(open(outDir+'final_consensus.fasta',"a")).write((open(outDir+sampleName+'.consensus.fasta',"r")).read())
		return False

def processClusterOutput(outDir,clusterFile, fastq, consensusFile, rankingMethod, pool, startIndex, endIndex, snpThresh, percentThresh, minorAlleleReadThresh, pvalThresh, maxSignatures, estimatedErrorRate, readCountSNPPvalThresh, readCountSNPVarFreqThresh, readCountSNPCoverageThresh, *optionalfilters):
	#numClusters=int((subprocess.check_output(['wc', '-l', clusterFile])).split(None, 1)[0])
	#09112017 startIndex=1
	#09112017 endIndex=getConsensusLength(consensusFile)
	with open(clusterFile) as f:
    		numClusters=sum(1 for _ in f)
	print str(numClusters) + " clusters"
	if numClusters > 1:
		with open(clusterFile,"r") as clusters: #where cluster output is the output of the cluster function
			print "Reading "+ clusterFile
			threads = []
			for line in clusters:
				print "Reading clusters"
				currCluster = line.strip()
				if os.path.isfile(currCluster):
					print currCluster
					currClusterFastq=currCluster+".fastq"
					getClusterFastq(fastq,currCluster,currClusterFastq)

					print "*****Recursive call for "+currCluster+"*****"
					#pool.spawn(pipeline,currClusterFastq,consensusFile,currCluster,rankingMethod, pool)
					pipeline(outDir,currClusterFastq,consensusFile,currCluster,rankingMethod, pool, startIndex, endIndex, snpThresh, percentThresh, minorAlleleReadThresh, pvalThresh, maxSignatures, estimatedErrorRate, readCountSNPPvalThresh, readCountSNPVarFreqThresh, readCountSNPCoverageThresh)
			clusters.close()
	else:
		print "Creating consensus for "+consensusFile+ "(processClusterOutput)"
		#command = 'cat <(echo ">'+consensusFile+'") <(tail -n\+2 '+consensusFile+'.consensus.fasta) >> outDir+final_consensus.fasta'
		#subprocess.call(command, shell = True,executable="/bin/bash")
		writeFinalClusterConsensus(outDir,consensusFile)
		#(open(outDir+'final_consensus.fasta',"a")).write((open(outDir+consensusFile+'.consensus.fasta',"r")).read())


def pipeline(outDir, initialFastq, referenceFile, sampleName, rankingMethod, pool, startIndex, endIndex, snpThresh, percentThresh,minorAlleleReadThresh, pvalThresh, maxSignatures, estimatedErrorRate, readCountSNPPvalThresh, readCountSNPVarFreqThresh, readCountSNPCoverageThresh, *optionalfilters):

	print ("*****Running pipeline for "+sampleName+"*****")

	start = time.time()

	try:
		os.makedirs(outDir)
	except OSError as exc:  # Python >2.5
		if exc.errno == errno.EEXIST and os.path.isdir(outDir):
			pass
		else:
			raise

	sampleNameSplit=os.path.split(sampleName)
	if len(sampleNameSplit)>1: sampleName=sampleNameSplit[1]

	if startIndex==0 or startIndex==None:
		startIndex=1
		endIndex=getConsensusLength(consensusFile)
	else:
		startIndex=int(startIndex)
		endIndex=int(endIndex)

	if snpThresh==0 or snpThresh==None:
		snpThresh=10
	else:
		snpThresh=int(snpThresh)



	print ("*****startIndex="+str(startIndex)+"*****")
	print ("*****endIndex="+str(endIndex)+"*****")
	print ("*****snpThresh="+str(snpThresh)+"*****")
	print ("*****percentThresh="+str(percentThresh)+"*****")


	# Declare fileNames
	initialBamfile =   outDir+str(sampleName) + ".bam"
	readCountsFile =  outDir+str(sampleName) + ".readcounts"
	consensusFile = outDir+str(sampleName) + ".consensus.fasta"
	realignedSampleName=str(sampleName) + ".realign"
	realignedBamfile =   outDir+str(realignedSampleName) + ".bam"
	realignedPileupFile =   outDir+str(realignedSampleName) + ".pileup"
	realignedSnpFile = outDir+str(realignedSampleName) + ".readcounts.snp"
	rankedFile = outDir+str(realignedSampleName) + ".rankedsnp"
	ranking_method = str(rankingMethod)
	clusterFile = outDir+str(sampleName) + "_clusters"

	# pvalThresh = 0.0001 #rankedSNPs  (default 10^-4)

	saved_args = locals()
    	print("pipeline parameters:", saved_args)


	consenseAndRealignSparc(initialFastq,referenceFile,initialBamfile, readCountsFile, consensusFile, realignedBamfile, estimatedErrorRate, readCountSNPPvalThresh, readCountSNPVarFreqThresh, readCountSNPCoverageThresh )

	rankSNVs(realignedSnpFile, rankedFile, rankingMethod, consensusFile, startIndex, endIndex, percentThresh, pvalThresh, minorAlleleReadThresh, maxSignatures)

	if cluster(outDir,realignedBamfile, rankedFile, sampleName, snpThresh):
		print("processClusterOutput next")
		processClusterOutput(outDir,clusterFile, initialFastq, consensusFile, rankingMethod, pool, startIndex, endIndex, snpThresh, percentThresh, minorAlleleReadThresh, pvalThresh, maxSignatures, estimatedErrorRate, readCountSNPPvalThresh, readCountSNPVarFreqThresh, readCountSNPCoverageThresh)

	end = time.time()
	print "Pipeline round for "+sampleName+" took " + str((end - start))+ "seconds"

def startPipeline(outDir,initialFastq, referenceFile, sampleName, rankingMethod, pool, startIndex, endIndex, snpThresh, percentThresh, minorAlleleReadThresh, pvalThresh, maxSignatures, estimatedErrorRate, readCountSNPPvalThresh, readCountSNPVarFreqThresh, readCountSNPCoverageThresh,  *optionalfilters):

	buffer=30
	print ("*****startPipeline for "+sampleName+"*****")
	referenceLength=getConsensusLength(referenceFile)
	if startIndex==None or int(startIndex)<1:
		startIndex=1
		print "*****Setting startIndex to 1"
	else:
		startIndex=int(startIndex)

	if endIndex==None or int(endIndex)<1 or int(endIndex)>referenceLength:
		endIndex=referenceLength
		print "*****Setting endIndex to "+str(referenceLength)
	else:
		endIndex=int(endIndex)

	start = time.time()
	m5File=outDir+sampleName+".initial.m5"
	firstRoundConsensusFile=os.path.splitext(referenceFile)[0]+'.firstRoundSparc.consensus.fasta'

	fullLengthReadIds=outDir+sampleName+"_fullLength"
	fullLengthReadFastq=outDir+sampleName+"_fullLength.fastq"
	nonFullLengthReadIds=outDir+sampleName+"_nonFullLength"
	#nonFullLengthReadFastq=outDir+sampleName+"_nonFullLength.fastq"

	saved_args = locals()
    	print("startPipeline parameters:", saved_args)


	alignFastqBlasrM5(initialFastq,referenceFile,m5File)# Align subreads to the initial reference and index resulting bam

	getFullLengthReads(outDir,sampleName,m5File,int(startIndex),int(endIndex),buffer)


	truncateFasta(referenceFile,referenceFile,int(startIndex),int(endIndex))
	getOrientedFastq(initialFastq,fullLengthReadIds,fullLengthReadFastq)
	#getOrientedFastq(initialFastq,nonFullLengthReadIds,nonFullLengthReadFastq)
	m5File=outDir+sampleName+".m5"
	alignFastqBlasrM5(fullLengthReadFastq,referenceFile,m5File)# Align subreads to the initial reference and index resulting bam
	getConsensusSparc(m5File,referenceFile,firstRoundConsensusFile)# Get consensus using Sparc (first round)

	#normalizedStartEnd=normalizeStartEnd(referenceFile,firstRoundConsensusFile,startIndex,endIndex)
	#startIndex=normalizedStartEnd[0]
	#endIndex=normalizedStartEnd[1]

	copyfile(fullLengthReadIds, outDir+sampleName)

	referenceLength=getConsensusLength(firstRoundConsensusFile)
	startIndex=buffer
	endIndex=referenceLength-buffer
	pipeline(outDir, fullLengthReadFastq, firstRoundConsensusFile, sampleName, rankingMethod, pool,startIndex,endIndex, snpThresh, percentThresh, minorAlleleReadThresh, pvalThresh, maxSignatures, estimatedErrorRate, readCountSNPPvalThresh, readCountSNPVarFreqThresh, readCountSNPCoverageThresh)
	writeFinalConsensusStats(outDir)
	end = time.time()
	print "Initialization took " + str((end - start))+ "seconds"

def main(argv):
   start = time.time()
   global writeUnaligned
   writeUnaligned = True
   outDir='/data/quasiSeqOut/multi/'
   initialFastq='/data/lu10mix_FL_cc.fastq'
   # referenceFile='/data/quasiSeqOut/201566f9-a5ec-4220-bde3-af488966e2b0/flu1PB.fa'
   referenceFile='/data/flu1PB.fa'
   sampleName='lu10mix_FL_cc'
   rankingMethod='percentage'
   pool=mp.Pool(16)
   startIndex=25
   endIndex=2270
   snpThresh='4'
   percentThresh=1.0
   minorAlleleReadThresh=10
   pvalThresh='0.01'
   maxSignatures='13000'
   estimatedErrorRate='0.01'
   readCountSNPPvalThresh=0.05
   readCountSNPVarFreqThresh=0.001
   readCountSNPCoverageThresh=10

   startPipeline(outDir,initialFastq, referenceFile, sampleName, rankingMethod, pool, startIndex, endIndex, snpThresh, percentThresh, minorAlleleReadThresh, pvalThresh, maxSignatures, estimatedErrorRate, readCountSNPPvalThresh, readCountSNPVarFreqThresh, readCountSNPCoverageThresh)
   pool.close()
   pool.join()
   #  for i in range(10):
   #      pool.apply_async(foo_pool, args = (i, ), callback = log_result)
   #  pool.close()
   #  pool.join()
   #  print(result_list)
   end = time.time()
   print "**********Full processing took " + str((end - start))+ "seconds**********"

if __name__ == "__main__":
   main(sys.argv[1:])

writeUnaligned = True
