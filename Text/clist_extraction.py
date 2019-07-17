###############################################################################
#Uses a short list of entries to extract from a much larger list of entries
#
#@author:charles.hefer@gmail.com
#@version:0.1
###############################################################################
import optparse
import csv


def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-l", "--list", default=None, dest="list",
					  help="The input file, short list")
	parser.add_option("-c", "--content", default=None, dest="content",
					  help="The content file, longer list")
	parser.add_option("-o", "--output", default=None, dest="output",
						help="The output file")
	parser.add_option("","--sep", default=" ", dest="sep",
						help="Specify the column seperator, default is ' '")
	parser.add_option("", "--unique-id", default="1", dest="col",
						help="The column to treat as the unique id, default is the first col (1)")
	(options, args) = parser.parse_args()
	
	if not options.list:
		parser.error("Need to specify the input list file")
	if not options.content:
		parser.error("Need to specify the input content file")
	try:
		column = int(options.col) - 1
	except:
		parser.error("Can not convert the --unique-id to an integer, try again")
		

	entries = []
	outlines = []
	unique_id = None


	with open(options.list, "r") as shortfile:
		for line in shortfile:
			line = line.rstrip()
			print(line)
			with open(options.content, "r") as longfile:
				for lline in longfile:
					cols = lline.split(options.sep)
					if cols[column] == line:
						outlines.append("\t".join([line, lline]))
						break
	
	with open(options.output, "w") as handle:
		for line in outlines:
			handle.write(line)

			
		
if __name__ == "__main__":
	__main__()
	