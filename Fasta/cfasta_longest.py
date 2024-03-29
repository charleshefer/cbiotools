###############################################################################
#cfasta_longest
#Reads in a fasta formatted file, writes out a the longest entry
#biopython
#
#@requires biopython
#@author:charles.hefer@gmail.com
#@version:0.1
###############################################################################
import optparse
from Bio import SeqIO


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
	
	#Read into RAM, will bomb if there are Gbs of entries
	#but this has not happened yet.
	#Better than IO overhead for most cases
	#fasta_records = []
	longest_record = {"record" : None,
					 "record_len" : 0}


	with open(options.input, "r") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			if len(record.seq) > longest_record["record_len"]:
				longest_record["record"] = record
				longest_record["record_len"] = len(record.seq)
				
	with open(options.output, "w") as outhandle:
		SeqIO.write([longest_record["record"]], outhandle, "fasta")
		
if __name__ == "__main__":
	__main__()