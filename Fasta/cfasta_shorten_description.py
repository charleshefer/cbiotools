###############################################################################
#cfasta_cleaner
#Reads in a fasta formatted file that contains several description seperated by
# a ">" (This is default if you extract from a blast databas
# Seperate the descriptions based on a keyword, and get rid of any duplicates
#biopython
#
#@requires biopython
#@author:charles.hefer@gmail.com
#@version:0.1
###############################################################################
import optparse, sys
from Bio import SeqIO


def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input file")
	parser.add_option("-o", "--ouput", default=None, dest="output",
					  help="The output file")
	parser.add_option("-k", "--keyword", default=None, dest="keyword",
					help="The keyword to search for")
	(options, args) = parser.parse_args()
	
	if not options.input:
		parser.error("Need to specify the input file")
	if not options.output:
		parser.error("Need to specify the output file")
	if not options.keyword:
		parser.error("Need to specify a keyword")
	
	#Read into RAM, will bomb if there are Gbs of entries
	#but this has not happened yet.
	#Better than IO overhead for most cases
	fasta_records = []
	dedub_db = []
	
	with open(options.input, "r") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			descriptions = record.description.split(">")
			short_description = None
			for description in descriptions:
				if options.keyword in description:
					short_description = description
					break
			if short_description: record.description = short_description
			else: record.description = descriptions[0]
			
			record.name = record.description.split(" ")[0]
				
			if record.name not in dedub_db:
				fasta_records.append(record)
				dedub_db.append(record.name)

	with open(options.output, "w") as outhandle:
		SeqIO.write(fasta_records, outhandle, "fasta")
		
if __name__ == "__main__":
	__main__()
	
