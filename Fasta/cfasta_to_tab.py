###############################################################################
#cfasta_to_tab
#Reads in a fasta formatted file, writes out a tab delimited file
#biopython
#
#@requires biopython
#@author:charles.hefer@gmail.com
#@version:0.1
###############################################################################
import optparse
from Bio import SeqIO


def parse_description(description, sep):
	"""returns a single text entry seperated by tabs"""
	return "\t".join(description.split(sep))

def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input file")
	parser.add_option("-s", "--sep", default="|", dest="sep",
					help="The seperator field for the description")
	parser.add_option("-o", "--ouput", default=None, dest="output",
					  help="The output file")
	(options, args) = parser.parse_args()
	
	if not options.input:
		parser.error("Need to specify the input file")
	if not options.output:
		parser.error("Need to specify the output file")
	
	#Read into RAM, will bomb if there are Gbs of entries
	#but this has not happened yet.
	#Better than IO overhead for most cases
	outlines = []
	
	with open(options.input, "r") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			outline = "%s\t%s\t%s" % (record.id, parse_description(record.description, options.sep), 
									str(record.seq))
			outlines.append(outline)
	
	with open(options.output, "w") as outhandle:
		for outline in outlines:
			outhandle.write(outline + "\n")
		
if __name__ == "__main__":
	__main__()
	