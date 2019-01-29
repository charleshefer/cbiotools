###############################################################################
#Convert genbank to fasta format
#
#@requires Biopython
#@author:charles.hefer@gmail.com
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
					  help="The output file")
	(options, args) = parser.parse_args()
	
	if not options.input:
		parser.error("Need to specify the input genbank file")
	if not options.output:
		parser.error("Need to specify the output fasta file")
		
	with open(options.input, "r") as handle:
		entries = [entry for entry in SeqIO.parse(handle, "genbank")]
	
	with open(options.output, "w") as handle:
		SeqIO.write(entries, handle, "fasta")
		
	
		
if __name__ == "__main__":
	__main__()
	