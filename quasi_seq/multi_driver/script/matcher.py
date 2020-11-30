from __future__ import division
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument("leaves" , help = "leaves file needed" )
parser.add_argument("fasta", help = "fasta file needed")
parser.add_argument("-of", "--outfastafile", default = "default", type = str)
parser.add_argument("-ol", "--outleavesfile", default = "default", type = str)

args = parser.parse_args()

leaves_file = str(args.leaves)
fasta_file = str(args.fasta)
outf_file = str(args.outfastafile)
outl_file = str(args.outleavesfile)

if outf_file == "default": 
	outf_file = "final_" + fasta_file
if outl_file == "default": 
	outl_file = "final_" + leaves_file 




finalf_file  = open(outf_file, "w+")
finalL_file = open(outl_file, "w+")

leaves_list = [] 
with open(leaves_file, "r") as leaves: 
	total_lines = 0

	for line in leaves:
		line = line.strip()  
		total_lines += sum(1 for line in open(line))
#finalL_file.write( str(total_lines)+ "\t"+ "TOTAL" + "\n")
with open(leaves_file, "r") as leave_file:
	for line  in leave_file:
		line= line.strip() 
	#	print "TEST" 
		num_lines = sum(1 for line in open(line))
		freq = num_lines/total_lines 
		freq = round(freq, 4)
		line_string = str(line) +".freq_" + str(freq)
	#	print "LINE STRING" +  line_string 
		leaves_list.append(line_string)
		finalL_file.write(str(num_lines) + "\t" + str(line) + "\n")
finalL_file.write(str(total_lines) + "\t" + "TOTAL")
print "total strains: " + str(len(leaves_list))		
print "THIS IS TOTAL LINES " + str(total_lines) 

with open(fasta_file, "r") as fasta: 
	counter = 0 
	for line in fasta: 
		if ">" in line:
			#`:wqprint "COUNTER " + str(counter)
			#print "COUNTER" + str(counter)
			#print  leaves_list
			line  = leaves_list[counter]
			finalf_file.write(">"+line + "\n")
			counter += 1
		else: 
			finalf_file.write(line) 
