###############################################################################
#cfasta_split
#Reads in a fasta formatted file, writes out chunks of the file
#biopython
#
#@requires biopython
#@author:charles.hefer@gmail.com
#@version:0.1
###############################################################################
import optparse, os
from Bio import SeqIO


def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input file")
	parser.add_option("-o", "--outdir", default=None, dest="outdir",
					  help="The output directory")
	parser.add_option("-c", "--chunks", default=10, dest="chunks",
						help="The number of files to split into")
	parser.add_option("-p", "--prefix", default="part", dest="prefix",
						help = "The prefix to the filename")
	(options, args) = parser.parse_args()
	
	if not options.input:
		parser.error("Need to specify the input file")
	if not options.outdir:
		parser.error("Need to specify the output directory")
	
	#Read into RAM, will bomb if there are Gbs of entries
	#but this has not happened yet.
	#Better than IO overhead for most cases
	fasta_records = []

	if not os.path.exists(options.outdir):
		os.mkdir(options.outdir)

	
	with open(options.input, "r") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			fasta_records.append(record)
	
	entries_per_file = len(fasta_records)/int(options.chunks)
	print(entries_per_file)
	file_counter = 1
	with open(options.outdir + options.prefix + "_" + str(file_counter) + ".fasta", "w") as outhandle:
		current_chunk = []
		for record in fasta_records:
			if len(current_chunk) < entries_per_file:
				current_chunk.append(record)
			else:
				SeqIO.write(current_chunk, outhandle, "fasta")
				current_chunk = [record]
				outhandle.close()
				file_counter = file_counter + 1
				outhandle = open(options.outdir + options.prefix + "_" + str(file_counter) + ".fasta", "w")
		SeqIO.write(current_chunk, outhandle, "fasta")
		outhandle.close()


if __name__ == "__main__":
	__main__()
	