###############################################################################
#cfasta_cleaner
#Reads in a fastq formatted file, writes out a corresponding pair file with
#NNs as the sequences, and zero low quality scores (!)
#biopython
#
#@requires biopython
#@author:charles.hefer@gmail.com
#@version:0.1
###############################################################################
import optparse
import gzip
from Bio import SeqIO
from Bio.Seq import Seq


def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input file")
	parser.add_option("-o", "--ouput", default=None, dest="output",
					  help="The output file")
	
	(options, args) = parser.parse_args()
	
	if not options.input:
		parser.error("Need to specify the input fastq file")
	if not options.output:
		parser.error("Need to specify the output fastq file")
	
	
	fastq_records = []
	with open(options.output, "w") as handle:
		SeqIO.write(fastq_records, handle, "fastq")

	
	with gzip.open(options.input, "rt") as handle:
		for record in SeqIO.parse(handle, "fastq"):
			record_description = record.description.split()
			record_description[-1] = record_description[-1].replace("1:", "2:")
			record.description = " ".join(record_description)

			seq = ["N" for base in record.seq]
			record.seq = Seq("".join(seq))
			
			quality = [float(0) for q in record.letter_annotations["phred_quality"]]
			record.letter_annotations["phred_quality"] = quality
			
			fastq_records.append(record)

			if len(fastq_records) == 100000:
				with open(options.output, "a") as handle:
					SeqIO.write(fastq_records, handle, "fastq")
				fastq_records = []

	with open(options.output, "a") as handle:
		SeqIO.write(fastq_records, handle, "fastq" )
	
		
if __name__ == "__main__":
	__main__()
	