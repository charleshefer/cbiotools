###############################################################################
#Download a fasta entry from the NCBI using EUTILS
#@requires:biopython
#@author:charles.hefer@agresearch.co.nz
#@version:0.1
###############################################################################
import optparse
from Bio import Entrez, SeqIO

Entrez.email="charles.hefer@gmail.com"

def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-a", "--accession", default=None, dest="accession",
					  help="The accession to retrieve")
	parser.add_option("-d", "--database", default=None, dest="database",
					help = "The database to retrieve from")
	parser.add_option("-o", "--output", default=None, dest="output",
					  help="The output file")
	(options, args) = parser.parse_args()
	
	if not options.accession:
		parser.error("Need to specify the input file")
	if not options.database:
		parser.error("Need to specify the database")
	if not options.output:
		parser.error("Need to specify the output file")
		
	handle = Entrez.efetch(db=options.database,
		id = options.accession,
		rettype="fasta")
		
	records = handle.read()
	#This is not a fasta object....

	with open(options.output, "w") as outhandle:
		outhandle.write(records)

if __name__ == "__main__":
	__main__()
	