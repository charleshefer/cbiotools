###############################################################################
#Extract partial matches from a fasta filename
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
					  help="The input file")
	parser.add_option("-o", "--ouput", default=None, dest="output",
					  help="The output file")
	parser.add_option("-m", "--match", default=None, dest="match",
					  help="The partial match for the fasta entry, case sensitive")
	parser.add_option("-f", "--file", default=None, dest="file",
					help = "A filename with a list of entries to match")
	(options, args) = parser.parse_args()
	
	
	if not options.input:
		parser.error("Need to specify the input file")
	if not options.output:
		parser.error("Need to specify the output file")
	if not options.match and not options.file:
		parser.error("Need to specify either some part of an entry name (-m) or a filename of entries (-f) to match")

	if options.file:
		with open(options.file, "r") as handle:
			entries = [e.rstrip() for e in handle]



	fasta_records = []
	with open(options.input, "r") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			if options.match:
				if options.match in record.name:
					fasta_records.append(record)
			if options.file:
				for entry in entries:
					if entry in record.name:
						fasta_records.append(record)
						break
				
	with open(options.output, "w") as outhandle:
		SeqIO.write(fasta_records, outhandle, "fasta")
	
if __name__ == "__main__":
	__main__()
	