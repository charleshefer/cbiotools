###############################################################################
#crockhopper_reformat
#Reads in a tab seperated rockhopper file (typically _transcripts.txt)
#and reformats the "predicted RNA field to a unique id"
#
#@author:charles.hefer@gmail.com
#@version:0.1
###############################################################################
import optparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input rockhopper file")
	parser.add_option("-o", "--ouput", default=None, dest="output",
					  help="The output fasta file")
	(options, args) = parser.parse_args()
	
	if not options.input:
		parser.error("Need to specify the input file")
	if not options.output:
		parser.error("Need to specify the output fasta file")
	
	records = []
	counter = 1
	
	with open(options.input, "r") as handle:
		for line in handle:
			line = line.rstrip()
			if "Sequence" in line: continue
			cols = line.split("\t")
			seq = Seq(cols[0])
			record = SeqRecord(seq, 
								id = "mRNA_%i" % counter,
								description = "")
			counter = counter + 1
			records.append(record)


	with open(options.output, "w") as outhandle:
		SeqIO.write(records, outhandle, "fasta")

if __name__ == "__main__":
	__main__()