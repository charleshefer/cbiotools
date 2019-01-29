###############################################################################
#Converts a ERCC fasta file to a ERCC gtf file
#
#@requires: biopython
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
					  help="The output gtf file")
	(options, args) = parser.parse_args()
	
	
	if not options.input:
		parser.error("Need to specify the input file")
	if not options.output:
		parser.error("Need to specify the output file")
		
	gtf_lines = []
	
	with open(options.input, "r") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			gtf_line = '%s\tercc\tgene\t1\t%s\t.\t+\t.\tgene_id "%s"; gene_version "1"; gene_name "%s"; gene_source "ercc"; gene_biotype "ercc" \n' % (record.name,
																		   len(record.seq) + 1,
																		   record.name,
																		   record.name)
			gtf_lines.append(gtf_line)
				
	with open(options.output, "w") as outhandle:
		for entry in gtf_lines:
			outhandle.write(entry)
	
if __name__ == "__main__":
	__main__()
	