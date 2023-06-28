###############################################################################
#Given an input list of samples, randomly assign the samples in a balanced
#way to different batches
#
# 
# @author:charles.hefer@agresearch.co.nz
#@version:0.1
###############################################################################
import optparse,pandas



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

	with open(options.input, "r") as handle:
		for line in handle:
			line = line.rstrip()
			if line.startswith("No.,"): #this is the column names
				#line = handle.next()
				number, sample_code, species, batch, processing, gastric, intestinal = line.split(",")
				print(number) 


if __name__ == "__main__":
	__main__()
	