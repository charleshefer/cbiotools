###############################################################################
#Converts a gtf file to bed file, while making the name column more friendly
#
#@author:charles.hefer@gmail.com
#@version:0.1
###############################################################################
import optparse
import sys


def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input gtf file")
	parser.add_option("-o", "--ouput", default=None, dest="output",
					  help="The output bed file")
	(options, args) = parser.parse_args()
	
	
	if not options.input:
		parser.error("Need to specify the input gtf file")
	if not options.output:
		parser.error("Need to specify the output bed file")
		
	bed_lines = []
	
	with open(options.input, "r") as handle:
		for line in handle:
			if line.startswith("#"): continue
			cols = line.split("\t")
			#print(cols)
			#print(cols[8].replace(" ", "="))
			#sys.exit()
			chromosome = cols[0]
			start, stop = cols[3:5]
			start = str(int(start) - 1)
			name = cols[8].replace(" ","=").replace(";=", "")
			bed_lines.append("\t".join([chromosome, start, stop, name]))
		
				
	with open(options.output, "w") as outhandle:
		outhandle.write("#BED format\n")
		for entry in bed_lines:
			outhandle.write(entry)
	
if __name__ == "__main__":
	__main__()
	