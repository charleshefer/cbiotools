###############################################################################
#cfasta_cleaner
#Reads in a fasta formatted file, writes out a standard format using
#biopython
#
#@requires biopython
#@author:charles.hefer@gmail.com
#@version:0.1
###############################################################################
import optparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
	fasta_records = []
	
	with open(options.input, "r") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			new_record = SeqRecord(record._seq[::-1], id=record.id, description="")
			fasta_records.append(new_record)
	
	with open(options.output, "w") as outhandle:
		SeqIO.write(fasta_records, outhandle, "fasta")
		
if __name__ == "__main__":
	__main__()
	