###############################################################################
#Amendst the name of a fasta sequence (appends to it)
#
#@requires: biopython
#@author:charles.hefer@gmail.com
#@version:0.1
###############################################################################
import optparse
from Bio import SeqIO
import sys

def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input file")
	parser.add_option("-o", "--ouput", default=None, dest="output",
					  help="The output file")
	parser.add_option("-s", "--string", default=None, dest="string",
					  help="Text to add to a fasta entry name")
	(options, args) = parser.parse_args()
	
	
	if not options.input:
		parser.error("Need to specify the input file")
	if not options.output:
		parser.error("Need to specify the output file")
	
	
	fasta_records = []
	with open(options.input, "r") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			record.id = record.id + options.string
			#print(record.id)
			record.id = record.id.replace('.p','')
			fasta_records.append(record)
			#break

	with open(options.output, "w") as outhandle:
		SeqIO.write(fasta_records, outhandle, "fasta")
	
if __name__ == "__main__":
	__main__()
	
