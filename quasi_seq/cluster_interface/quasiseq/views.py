# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.http import HttpResponse
from django.shortcuts import render,redirect
from random import randint
import subprocess
from collections import defaultdict
#import gevent.monkey
#gevent.monkey.patch_all()
#import gevent
from multiprocessing import Pool
import multiprocessing.pool
import errno
from celery.result import AsyncResult
from django.conf import settings

import calendar
import time
import os

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
#from driver import fullrecursionpipeline
#from sparc import pipeline
#from linearQuasiSeq import pipeline
from tasks import pipeline_driver
#from driver import mediary

#from linear_sparc import fullrecursionpipeline, mediary
from multiprocessing import Pool, freeze_support, Manager

# Create your views here.
def index(request):
	#return render(request, "quasiseq/chart.html")
	response = render(request, "quasiseq/intro_form.html")

	var = subprocess.Popen("nproc", stdout = subprocess.PIPE)
	nproc = var.stdout.read() # number of cores on machine
#	pool = gevent.pool.Pool(int(nproc))

	#test(session_id)

	#currentTask=None
	if(request.method == "POST"):

		if(request.POST.get("Pro4mix")):
			return redirect("Pro4mix")



		session_id="quasiSeq_"+str(int(calendar.timegm(time.gmtime())))
		outDir="/data/quasiSeqOut/" #+session_id+"/"
		tmpDir="/data/quasiSeqOut/tmp/"
		try:
			os.makedirs(outDir)
		except OSError as exc:  # Python >2.5
			if exc.errno == errno.EEXIST and os.path.isdir(outDir):
				pass
			else:
				raise
		try:
			os.makedirs(tmpDir)
		except OSError as exc:  # Python >2.5
			if exc.errno == errno.EEXIST and os.path.isdir(tmpDir):
				pass
			else:
				raise

		print ("Session key:"+session_id)
		#num = randint(100000,999999)
		num = 99999
		#sess = request.session.session_key
		#request.session["ID"] = num
		uploaded_file = request.FILES['my_file']
		sampleName = request.FILES["my_file"].name
		sampleName = sampleName.strip(".fastq")
		samplefilename = sampleName + ".fastq"

		with open(tmpDir+samplefilename , "wb+") as samplefile: # add the plus so that it will create the file if it does not exist
			# get the session ID and append that to sampleID.txt file name so that it is a different file for every user and the same file does not get overwritten
			for chunk in uploaded_file.chunks():
				samplefile.write(chunk)

		reference_file = request.FILES["reference_file"]
		uploaded_reference = request.FILES["reference_file"].name
		#uploaded_reference = uploaded_reference.strip("fasta")
		referencefilename= uploaded_reference

		with open(tmpDir+referencefilename, "wb+") as  referencefile :
			for chunk in reference_file.chunks():
				referencefile.write(chunk)

		rankingMethod = request.POST.get("rankingMethod")
		if rankingMethod == None:
			rankingMethod = 'percentage'


		#fullrecursionpipeline(filename, reference, ranking, startIndex, endIndex, percentThresh, pvalThresh, minorAlleleReadThresh, snpThresh, False, pool)
		startIndex = request.POST["start"]
		endIndex = request.POST["end"]

		if startIndex == None or startIndex == '':
			startIndex = 0

		if endIndex == None or endIndex == '':
			endIndex = None

		percentThresh = request.POST["percentThresh"]
		pvalThresh = request.POST["pvalThresh"]
		minorAlleleReadThresh = request.POST["minorAlleleReadThresh"]
		snpThresh = request.POST.get("snpThresh")
		estimatedErrorRate = request.POST.get("estimatedErrorRate")
		maxSignatures = request.POST.get("maxSignatures")

		# Default threshold for ReadCount_snp.java
		readCountSNPPvalThresh=0.05
		readCountSNPVarFreqThresh=0.001
		readCountSNPCoverageThresh=10


		if percentThresh == None:
			percentThresh = 3.0
		if pvalThresh == None:
			pvalThresh = 0.001
		if minorAlleleReadThresh == None:
			minorAlleleReadThresh = 30
		if snpThresh == None:
			snpThresh = 40
		if estimatedErrorRate == None:
			estimatedErrorRate = 0.01
		if maxSignatures == None:
			maxSignatures = 12000

		saved_args = locals()
	    	print("index view parameters:", saved_args)

		currentTask=pipeline_driver.delay(outDir,sampleName, referencefilename, rankingMethod, startIndex, endIndex, percentThresh, pvalThresh, minorAlleleReadThresh, snpThresh, maxSignatures, estimatedErrorRate, readCountSNPPvalThresh, readCountSNPVarFreqThresh, readCountSNPCoverageThresh)

		return redirect("result", currentTask.task_id, sampleName)
	else:
		return response


def Pro4mix(request):
	response = render(request, "quasiseq/pro4_sub80.html")
	session_id="quasiSeq_"+str(int(calendar.timegm(time.gmtime())))
	outDir="/data/quasiSeqOut/" #+session_id+"/"
	tmpDir="/data/quasiSeqOut/tmp/"

	if (request.method == "POST"):

		if(request.POST.get("Pro4mix")):
			return redirect("Pro4mix")

		if (request.POST.get("home")):
			return redirect("index")
		try:
			os.makedirs(outDir)
		except OSError as exc:  # Python >2.5
			if exc.errno == errno.EEXIST and os.path.isdir(outDir):
				pass
			else:
				raise
		try:
			os.makedirs(tmpDir)
		except OSError as exc:  # Python >2.5
			if exc.errno == errno.EEXIST and os.path.isdir(tmpDir):
				pass
			else:
				raise

		print ("Session key:"+session_id)
		uploaded_file = settings.EXAMPLE_PATH+"/examples/pro4_sub80.fastq"
		sampleName = "pro4_sub80.fastq"
		sampleName = sampleName.strip(".fastq")
		samplefilename = sampleName + ".fastq"
		opened_file = open(uploaded_file, "r")
		with open(tmpDir+samplefilename , "wb+") as samplefile: # add the plus so that it will create the file if it does not exist
			# get the session ID and append that to sampleID.txt file name so that it is a different file for every user and the same file does not get overwritten
			for line in opened_file:
				samplefile.write(line)

		reference_file = settings.EXAMPLE_PATH+"/examples/HXB2.fasta"
		uploaded_reference = "HXB2.fasta"
		#uploaded_reference = uploaded_reference.strip("fasta")
		referencefilename= uploaded_reference
		opened_ref_file = open(reference_file, "r")

		with open(tmpDir+referencefilename, "wb+") as  referencefile :
			for line in opened_ref_file:
				referencefile.write(line)

		rankingMethod = request.POST.get("rankingMethod")
		if rankingMethod == None:
			rankingMethod = 'percentage'
		#fullrecursionpipeline(filename, reference, ranking, startIndex, endIndex, percentThresh, pvalThresh, minorAlleleReadThresh, snpThresh, False, pool)
		startIndex = request.POST.get("start")
		endIndex = request.POST.get("end")

		if startIndex == None or startIndex == '':
			startIndex = 0

		if endIndex == None or endIndex == '':
			endIndex = None

		percentThresh = request.POST["percentThresh"]
		pvalThresh = request.POST["pvalThresh"]
		minorAlleleReadThresh = request.POST["minorAlleleReadThresh"]
		snpThresh = request.POST.get("snpThresh")
		estimatedErrorRate = request.POST.get("estimatedErrorRate")
		maxSignatures = request.POST.get("maxSignatures")

		# Default threshold for ReadCount_snp.java
		readCountSNPPvalThresh=0.05
		readCountSNPVarFreqThresh=0.001
		readCountSNPCoverageThresh=10

		if percentThresh == None:
			percentThresh = 3.0
		if pvalThresh == None:
			pvalThresh = 0.001
		if minorAlleleReadThresh == None:
			minorAlleleReadThresh = 30
		if snpThresh == None:
			snpThresh = 40
		if estimatedErrorRate == None:
			estimatedErrorRate = 0.01
		if maxSignatures == None:
			maxSignatures = 12000

		saved_args = locals()
	    	print("index view parameters:", saved_args)

		currentTask=pipeline_driver.delay(outDir,sampleName, referencefilename, rankingMethod, startIndex, endIndex, percentThresh, pvalThresh, minorAlleleReadThresh, snpThresh, maxSignatures, estimatedErrorRate, readCountSNPPvalThresh, readCountSNPVarFreqThresh, readCountSNPCoverageThresh)

		return redirect("result", currentTask.task_id, sampleName)
	else:
		return response


def result(request,taskId, sampleName):

	taskResult=AsyncResult(taskId)
	taskState=taskResult.state
	statusbul = False
	specieslist = []
	domain = request.get_host()
	if (taskState == "SUCCESS") :
		statusbul = True
		"""
		path =  "/data/quasiSeqOut/" + str(taskId) + "/" + "final_consensus.fasta"

		filename = open(path, "r")
		response = HttpResponse(filename, content_type = "text/csv")
		response["Content-Disposition"] = "attachment; filename = final_consensus.fasta"

		return response
		"""
		specieslist = []
		labellist = []
		valuelist = []
		counter = 0
		with open("/data/quasiSeqOut/" + str(taskId) + "/" + "final_consensus_stats.txt", "r") as readfile :
			for line in readfile:
				linelist = line.split(",")
				specieslist.append(linelist)
				counter +=1

		for lister in specieslist:
				labellist.append(str(str((lister[1]).strip("\n")).strip("'")))
				valuelist.append(int(lister[0]))
		labellist = labellist[:-1] # here we are removing the last element that is just the total
 		valuelist = valuelist[:-1]

		#labellist = ["10clones1", "10clones2", "10clones3"]


		# if(request.method == "POST"):
		# 	return redirect("downloadFile", taskId, sampleName)

		response = render(request, "quasiseq/results.html", {"taskState" : taskState, "taskId" : taskId, "statusbul" : statusbul, "sampleName" : sampleName, "specieslist": specieslist, "labellist" : labellist,  "valuelist" : valuelist })
		return response

	url = request.get_full_path()
	return HttpResponse("<html><script>setTimeout(function(){window.location.reload(1);}, 10000);</script>\nTask status is: "+taskState+" <br><br>Stay in this page for your results. The page will refresh automatically every 10 seconds.<br><br>Or you can check for results later with this url: " + domain + url + " </html>" )


def downloadConsensus(request, taskId, sampleName):
		path =  "/data/quasiSeqOut/" + str(taskId) + "/" + "final_consensus.fasta"

		filename = open(path, "r")
		response = HttpResponse(filename, content_type = "text/csv")
		outFileName=str(taskId)+"_"+sampleName+"_final_consensus.fasta"
		response["Content-Disposition"] = "attachment; filename = "+outFileName

		return response

def downloadFrequency(request, taskId, sampleName):
		path =  "/data/quasiSeqOut/" + str(taskId) + "/" + "final_consensus_stats.txt"

		filename = open(path, "r")
		response = HttpResponse(filename, content_type = "text/csv")
		outFileName=str(taskId)+"_"+sampleName+"_final_consensus_stats.fasta"
		response["Content-Disposition"] = "attachment; filename = "+outFileName

		return response



	#	ref_file = request.FILES['reference_file']
