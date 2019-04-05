###############################################################################
#Parses the last column of a gtf file
#
#@author:charles.hefer@agresearch.co.nz
#@version:0.1
###############################################################################
import optparse
import re


def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input gtf file")
	parser.add_option("-o", "--ouput", default=None, dest="output",
					  help="The output csv delimited file")
	(options, args) = parser.parse_args()
	
	if not options.input:
		parser.error("Need to specify the input file")
	if not options.output:
		parser.error("Need to specify the output file")

	transcript_id_re = re.compile('transcript_id "(.+?)";')
	gene_id_re = re.compile('.+?gene_id "(.+?)";')
	gene_name_re = re.compile('.+?gene_name "(.+?)";')
	

	transcripts = dict()

	with open(options.input, "r") as handle:
		for line in handle:
			cols = line.split("\t")
			if len(cols) < 9: pass
			if transcript_id_re.match(cols[8]):
				transcript_id = transcript_id_re.match(cols[8]).groups()[0]
			else: transcript_id = "NA"
			if gene_id_re.match(cols[8]):
				gene_id = gene_id_re.match(cols[8]).groups()[0]
			else: gene_id = "NA"
			if gene_name_re.match(cols[8]):
				gene_name = gene_name_re.match(cols[8]).groups()[0]
			else: gene_name = "NA"
			if transcript_id not in transcripts.keys():
				transcripts[transcript_id] = {"transcript_id": transcript_id,
											"gene_id" : gene_id,
											"gene_name" : gene_name}
	with open(options.output, "w") as handle:
		handle.write("transcript_id,gene_id,gene_name\n")
		for entry in transcripts.keys():
			handle.write("%s,%s,%s\n" % (transcripts[entry]["transcript_id"],
											transcripts[entry]["gene_id"],
											transcripts[entry]["gene_name"]))		
		
if __name__ == "__main__":
	__main__()
	