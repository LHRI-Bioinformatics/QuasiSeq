
import os
import argparse
from collections import defaultdict
from operator import itemgetter

parser = argparse.ArgumentParser()
parser.add_argument("filename" , help = "input file needed" )
parser.add_argument("outputfile",help = "output file needed" )
parser.add_argument("ranking_method", choices = ["pval", "percentage"], default = "percentage", help = "first preference ranking method")
parser.add_argument("referencefile", help = "reference file needed") # take this line out if you are running this script by itself

# these next arguments are all optional filters . There are default values if any of these are not specified by the user
parser.add_argument("-s","--starting_index",type = int, default = 0, help = "optional starting index")
parser.add_argument("-e","--ending_index", type = int, default = None, help = "optional ending index")
parser.add_argument("-perc","--percentage_thresh", type = float, default = 3.0, help = "optional lowest threshold for percentage")
parser.add_argument("-pval","--pval_thresh", type = float, default = .001, help = "optional highest threshold for pvalue")
parser.add_argument("-r","--read2_thresh", type = int, default = 10, help = "optional lowest threshold for the second read")
parser.add_argument("-numSig","--number_of_signatures", type = int, default = 12000, help = "optional cap number Of signatures")
args = parser.parse_args()
snipDic = defaultdict(list)
listtuples = []
filename = args.filename
reference = args.referencefile 
startingindex = args.starting_index
percentagethresh = args.percentage_thresh
pvalthresh = args.pval_thresh
#pvalthresh = 0.01
read2thresh = args.read2_thresh
numSig = args.number_of_signatures

# ending index is defined later because there is no default value if it is non-existent 
with open(filename,"r") as inputfile:
    headingsline = inputfile.readline()
    headingsline = headingsline.strip(" ")
    for index, line in enumerate (inputfile):
    	#print line
        newline = line.split() # this is a list of all the elements in each line of the text file
        indexKey = index
        snipDic[index] = newline
rankmethod = args.ranking_method

for dicIndex,key in enumerate(snipDic):
    if (int(snipDic[key][1]) < startingindex):
        continue  # skip this if its out of bounds

    if args.ending_index:
        endingindex = int(args.ending_index)
        if (int(snipDic[key][1]) > endingindex):
            continue

    percentage = float(snipDic[key][6].strip('%'))
    pval = float(snipDic[key][7])
    read2 = int(snipDic[key][5])
#    if(percentage > percentagethresh):
#    	print percentage
#    if(pval < pvalthresh):
#    	print pval
#    if(read2 > read2thresh):
#    	print read2
    
    if(percentage < percentagethresh) or (pval > pvalthresh) :
        continue # skip this if its out of the threshold limits

    if(read2 < read2thresh):
        continue

    tup = (snipDic[key][1], percentage, pval, dicIndex)
   # print (tup)
    listtuples.append(tup)
	
if rankmethod == "percentage":
    listtuples = sorted(listtuples, key = itemgetter(2))
    sortedtuples = sorted(listtuples, key = itemgetter(1), reverse  = True)
  #  print (listtuples)
 #   print (sortedtuples)

else:
    listtuples = sorted(listtuples, key = itemgetter(1), reverse = True)
    sortedtuples = sorted(listtuples, key = itemgetter(2))
#print (sortedtuples)
i = 0

output_file = args.outputfile
f = open(output_file, 'w')
f.close()

outputfile = open(output_file,"a")
outputfile.write(str(len(sortedtuples)) + "\t" +str(numSig) + "\t(selected SNVs with reference "+ reference+ " and capNumberOfSignatures)\n")
#outputfile.write(headingsline)

for tup in sortedtuples:
 #    i +=1
#     print i
     lister = snipDic[tup[3]][1]
     result = "".join(lister)
     outputfile.write(result)
     outputfile.write("\n")
outputfile.close()
