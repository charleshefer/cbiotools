###############################################################################
#Creates a header files from the fasta entries in ENSEMBL format
#
#@requires: biopython
#@author:charles.hefer@gmail.com
#@version:0.1
###############################################################################
import optparse
from Bio import SeqIO
import sys
import re

def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input fasta file")
	parser.add_option("-o", "--output", default=None, dest="output",
					  help="The output file")
	(options, args) = parser.parse_args()
	
	header_re = re.compile(".+?[chromosome|scaffold]:(?P<position>.+?)\sgene:(?P<gene>.+?)\s.+?(gene_symbol:(?P<gene_symbol>.+?)\s)*(description:(?P<description>.+?))*$")
	
	if not options.input:
		parser.error("Need to specify the input file")
	if not options.output:
		parser.error("Need to specify the output file")
		
	#first, parse out the fasta name and description
	
	outlines = []
	with open(options.input, "r") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			outline = [record.id]
			
			re_match = header_re.match(record.description)
			if re_match:
				for match in ["position", "gene", "gene_symbol", "description"]:
					if re_match.group(match):
						outline.append(re_match.group(match))
					else:
						outline.append("NA")
			else:
				print(record.description)
				sys.exit()
			outlines.append("\t".join(outline))
	
	with open(options.output, "w") as outhandle:
		for line in outlines:
			outhandle.write(line + "\n")
	
if __name__ == "__main__":
	__main__()
	