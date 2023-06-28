###############################################################################
#EMBOSS:
#format the output from restrict into a bed  file
#
#parameters used:
#	restrict -sequence REFERENCE -enzymes BamHI,BglII,EcoRI,HindIII,PstI,SmaI -outfile test_restrict.out
#@author:charles.hefer@agresearch.co.nz
#@version:0.1
###############################################################################
import optparse
from Bio import Restriction
from Bio.Restriction import *
from Bio import SeqIO


def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input restrict file")
	parser.add_option("-e", "--enzymes", default="EcoRI,HindIII", dest='enzymes',
                      help="Comma separated list of enzymes")
	parser.add_option("-o", "--ouput", default=None, dest="output",
					  help="The output bed file")
	(options, args) = parser.parse_args()
	
	if not options.input:
		parser.error("Need to specify the input restrict file")
	if not options.output:
		parser.error("Need to specify the output bed file")
		
	#build the batch analysis
	rb = RestrictionBatch(options.enzymes.split(","))
	print("Searching for these enxymes: ")
	print(rb)
	print("")


	enzyme_lengths = {}
	for enzyme in rb:
		enzyme_lengths[enzyme] = len(enzyme.site)		
	print enzyme_lengths
   
	with open(options.input, "r") as handle:
		with open(options.output, "w") as outhandle:
			for record in SeqIO.parse(handle, "fasta"):
				result = Analysis(rb, record.seq, linear=False)
				for key in result.mapping.keys():
						#print(key)
						#print(result.mapping[key])
						for position in result.mapping[key]:
							outhandle.write("%s\t%s\t%s\t%s\n" % (record.name, 
                                int(position) - 1, 
                                enzyme_lengths[key] + int(position) - 1 - 1, 
                                str(key) + " cut site"))
				#print(record.name, )
   #print(result.full())
   





			


if __name__ == "__main__":
	__main__()
	