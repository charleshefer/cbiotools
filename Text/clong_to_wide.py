###############################################################################
#Takes a long format document, merges entries and return a wide format
#
#@author:charles.hefer@gmail.com
#@version:0.1
###############################################################################
import optparse
import csv


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
		column = int(options.col) - 1
	except:
		parser.error("Can not convert the --unique-id to an integer, try again")
		

	entries = []
	outlines = []
	unique_id = None


	with open(options.input, "r") as csv_file:
		handle = csv.reader(csv_file)
		for line in handle:
			
			#line = line.rstrip()
			#print(line)
			#cols = line.split(options.sep)
			cols = line
			#print(cols[0])
			if len(cols) == 0:
				#print("No columns in file, or there are empty rows. Exiting")
				sys.exit()
			
			

			if cols[column] not in entries and unique_id is None:
				#print("First time I see this entry:", cols[column])
				unique_id = cols[column]
				entries.append(unique_id)
				col_dict = {}
				col_index = 0
				for col in cols:
					col_dict[col_index] = [col]
					col_index = col_index + 1
				#print(col_dict)
				continue
			
			#print("Before Second")
			if (cols[column] in entries) and (unique_id == cols[column]):
				#continue adding
				#print("Inside Second")
				col_index = 0
				#print(col_dict)
				#print(len(cols))
				for col in cols:
					col_dict[col_index].append(col)
					col_index = col_index + 1
				#print(col_dict)
				continue

			#print("Before third")
			#print(cols[column], entries, unique_id)
			if (cols[column] not in entries) and (unique_id is not None):
				#print("Inside Third")
				outline = []
				for i in range(0, len(col_dict.keys())):
					col_dict[i] = "; ".join(list(set(col_dict[i])))
					outline.append(col_dict[i])
				outlines.append("\t".join(outline))	
								
				unique_id = cols[column]
				entries.append(unique_id)
				col_dict = {}
				col_index = 0
				for col in cols:
					col_dict[col_index] = [col]
					col_index = col_index + 1
				continue
			
	outline = []
	for i in range(0, len(col_dict.keys())):
		col_dict[i] = "; ".join(list(set(col_dict[i])))
		outline.append(col_dict[i])
	outlines.append("\t".join(outline))
	
	
	with open(options.output, "w") as handle:
		for line in outlines:
			handle.write(line + "\n")




			
		
if __name__ == "__main__":
	__main__()
	