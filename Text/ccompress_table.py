###############################################################################
#Takes a long format document, merges entries and return a wide format
#
#@author:charles.hefer@gmail.com
#@version:0.1
###############################################################################
import optparse
import csv
import copy


def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input file")
	parser.add_option("-o", "--ouput", default=None, dest="output",
					  help="The output file")
	parser.add_option("","--sep", default=",", dest="sep",
						help="Specify the column seperator, default is ','")
	parser.add_option("", "--unique-id", default="1", dest="col",
						help="The column to treat as the unique id, default is the first col (1)")
	(options, args) = parser.parse_args()
	
	if not options.input:
		parser.error("Need to specify the input file")
	if not options.output:
		parser.error("Need to specify the output file")
	try:
		ucol = int(options.col) - 1
	except:
		parser.error("Can not convert the --unique-id to an integer, try again")
	if options.sep == "\\t":
		options.sep = "\t"


	unique_lines = {}

	with open(options.input, "r") as handle:
		for line in handle:
			line = line.rstrip()
			cols = line.split(options.sep)
			print(cols)
			
			#if unqique column has neven been seen before:
			if cols[ucol] not in unique_lines.keys():
				unique_lines[cols[ucol]] = {}
				for i in range(0,len(cols)):
					unique_lines[cols[ucol]][i] = [cols[i]]
			else:
				for i in range(0,len(cols)):
					#this can through a keyerror if
					#the first instance of the ucol has no
					#other columns
					try:
						unique_lines[cols[ucol]][i].append(cols[i])
					except KeyError:
						unique_lines[cols[ucol]][i] = [cols[i]]
	for entry in unique_lines.keys():
		for column in unique_lines[entry]:
			unique_lines[entry][column] = list(set(unique_lines[entry][column]))
	
	#write to file
	with open(options.output, "w") as handle:
		for entry in unique_lines.keys():
			outline = []
			for column in unique_lines[entry]:
				outline.append(";".join(unique_lines[entry][column]))
			handle.write(options.sep.join(outline) + "\n")
	
if __name__ == "__main__":
	__main__()
	