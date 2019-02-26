###############################################################################
#Extract annotation from the description fields of a fasta file
#
#@requires: biopython
#@author:charles.hefer@gmail.com
#@version:0.1
###############################################################################
import optparse
from Bio import SeqIO
import sys

def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input fasta file")
	parser.add_option("-o", "--ouput", default=None, dest="output",
					  help="The output file")
	parser.add_option("-f", "--file", default=None, dest="file",
					help = "A filename with a list of entries to match")
	parser.add_option("-c", "--accession_col", default="1", dest="accession_col",
					help = "The column in the entries file to use as search ('1' based), default is the first column (1)")
	parser.add_option("-s", "--substitute", default=None, dest="substitute",
					help = "Substitute this character in the list of queries, default is None")
	parser.add_option("-w", "--sub_with", default=None, dest="substitute_with",
					help = "Substitute the -s character with this character, default is None")
	(options, args) = parser.parse_args()
	
	
	if not options.input:
		parser.error("Need to specify the input file")
	if not options.output:
		parser.error("Need to specify the output file")
	if not options.file:
		parser.error("A filename of entries (-f) to match")
	try:
		options.accession_col = int(options.accession_col) - 1
	except Exception:
		parser.error("The -c option needs to be an integer")

	
	#first, parse out the fasta name and description
	fasta_entries = {}
	with open(options.input, "r") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			record_name = "|".join(record.name.split("|")[:2])
			fasta_entries[record_name] = record.description.split("|")[-1]

	outlines = []
	if options.file:
		with open(options.file, "r") as handle:
			for line in handle:
				line = line.rstrip()
				cols = line.split("\t")
				if options.substitute and options.substitute_with:
					accession = cols[options.accession_col].replace(options.substitute, options.substitute_with)
				else:
					accession = cols[options.accession_col]
				try:
					cols.append(fasta_entries[accession])
				except KeyError:
					pass
				outlines.append("\t".join(cols))

	with open(options.output, "w") as outhandle:
		for line in outlines:
			outhandle.write(line + "\n")
	
if __name__ == "__main__":
	__main__()
	