###############################################################################
#split a fasta file into multiple entries
#
#@author:charles.hefer@agresearch.co.nz
#@version:0.1
###############################################################################
import optparse
from Bio import SeqIO


def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input file")
	parser.add_option("-o", "--ouput", default=None, dest="output",
					  help="The output directory")
	(options, args) = parser.parse_args()
	
	if not options.input:
		parser.error("Need to specify the input file")
	if not options.output:
		parser.error("Need to specify the output directory")
		
	with open(options.input,"r") as handle:
		for entry in SeqIO.parse(handle, "fasta"):
			with open(options.output + "/" + entry.id + ".fasta", "w") as outhandle:
				SeqIO.write(entry, outhandle, "fasta")

if __name__ == "__main__":
	__main__()
	