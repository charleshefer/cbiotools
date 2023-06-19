###############################################################################
#crockhopper_reformat
#Reads in a tab seperated rockhopper file (typically _transcripts.txt)
#and reformats the "predicted RNA field to a unique id"
#
#@author:charles.hefer@gmail.com
#@version:0.1
###############################################################################
import optparse


def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input rockhopper file")
	parser.add_option("-o", "--ouput", default=None, dest="output",
					  help="The output rockhopper file")
	(options, args) = parser.parse_args()
	
	if not options.input:
		parser.error("Need to specify the input file")
	if not options.output:
		parser.error("Need to specify the output file")
	
	records = []
	counter = 1
	
	with open(options.input, "r") as handle:
		for line in handle:
			line = line.rstrip()
			cols = line.split("\t")
			if len(cols) == 1: break
			if cols[6] == "predicted RNA":
				counter = counter + 1
				cols[6] = cols[6].replace(" ", "_") + "_" + str(counter)
			new_line = "\t".join(cols)
			records.append(new_line)


	with open(options.output, "w") as outhandle:
		for record in records:
			outhandle.write(record + "\n")
	
if __name__ == "__main__":
	__main__()