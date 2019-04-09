###############################################################################
#Report the length of NNs in the fasta file, by entry 
#
#@author:charles.hefer@agresearch.co.nz
#@version:0.1
###############################################################################
import optparse
from Bio import SeqIO

import re

def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input file")
	parser.add_option("-o", "--ouput", default=None, dest="output",
					  help="The output file")
	(options, args) = parser.parse_args()

	if not options.input:
		parser.error("Need to specify the input file")
	if not options.output:
		parser.error("Need to specify the output file")
	
	space_re = re.compile("[A,C,G,T]*(N+)[A,C,G,T]*")
	records = []


	with open(options.output, "w") as outhandle:
		with open(options.input, "r") as handle:
			for entry in SeqIO.parse(handle, "fasta"):
				res = space_re.findall(str(entry.seq))
				if res:
					print(entry.id)
					for r in res:
						outhandle.write("%s\t%s\n" % (entry.id, len(r)))
		
if __name__ == "__main__":
	__main__()
	