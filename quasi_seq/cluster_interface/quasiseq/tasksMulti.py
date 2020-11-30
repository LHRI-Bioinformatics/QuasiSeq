from celery.decorators import task
# from linearQuasiSeq import pipeline
#from linearQuasiSeq import alignFastqBlasrM5
#from linearQuasiSeq import getConsensusSparc
# from linearQuasiSeq import startPipeline
from multiQuasiSeq import startPipeline
import os


@task(name = "pipeline_driver")
def pipeline_driver(outDir,sampleName, referencefilename, rankingMethod, startIndex, endIndex, percentThresh, pvalThresh, minorAlleleReadThresh, snpThresh, maxSignatures, estimatedErrorRate, readCountSNPPvalThresh, readCountSNPVarFreqThresh, readCountSNPCoverageThresh):
	# minorAlleleReadThresh=read2thresh

	print("**************Calling Pipeline*************")
	print("************My task id is: "+pipeline_driver.request.id+"***************")

	outDir=outDir+pipeline_driver.request.id+"/"

	try:
		os.makedirs(outDir)
	except OSError as exc:  # Python >2.5
		if exc.errno == errno.EEXIST and os.path.isdir(outDir):
			pass
		else:
			raise

	os.rename("/data/quasiSeqOut/tmp/"+sampleName+".fastq", outDir+sampleName+".fastq")
	os.rename("/data/quasiSeqOut/tmp/"+referencefilename, outDir+referencefilename)

	if percentThresh==None or percentThresh==0:
		percentThresh=1
	else:
		percentThresh=float(percentThresh)

	if minorAlleleReadThresh==None or minorAlleleReadThresh==0:
		minorAlleleReadThresh=10
	else:
		minorAlleleReadThresh=int(minorAlleleReadThresh)

	referenceFile=outDir+referencefilename
	initialFastq=outDir+sampleName+".fastq"

	pool=[]

	saved_args = locals()
    	print("pipeline_driver parameters:", saved_args)

	startPipeline(outDir,initialFastq, referenceFile, sampleName, rankingMethod, pool, startIndex, endIndex, snpThresh, percentThresh, minorAlleleReadThresh, pvalThresh, maxSignatures, estimatedErrorRate, readCountSNPPvalThresh, readCountSNPVarFreqThresh, readCountSNPCoverageThresh)

	#pipeline(outDir,outDir+uploaded_name+".fastq", outDir+referencefilename, uploaded_name, ranking_method, None, startIndex, endIndex, snpThresh )
