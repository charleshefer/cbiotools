###############################################################################
#Jessica Gathercole
#Annotation of peptides
#
#@author:charles.hefer@agresearch.co.nz
#@version:0.1
###############################################################################
import optparse


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
		
if __name__ == "__main__"
	__main__()
	