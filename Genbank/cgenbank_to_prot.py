###############################################################################
#Convert genbank to protein sequence
#
#@requires Biopython
#@author:charles.hefer@gmail.com
#@version:0.1
###############################################################################
import optparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
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
		parser.error("Need to specify the input genbank file")
	if not options.output:
		parser.error("Need to specify the output fasta file")
	
	entries = []
	with open(options.input, "r") as handle:
		for seq_record in SeqIO.parse(handle, "genbank"):
			for seq_feature in seq_record.features:
				#print(seq_feature.type)
				if seq_feature.type == "CDS":
						assert len(seq_feature.qualifiers['translation']) == 1
						#print(seq_feature.qualifiers)
						prot_seq = SeqRecord(Seq(seq_feature.qualifiers['translation'][0]))
						prot_seq.id = seq_feature.qualifiers['locus_tag'][0]
						if "protein_id" in seq_feature.qualifiers.keys():
							prot_seq.description = "%s|%s|%s" % (seq_record.name, seq_feature.qualifiers["protein_id"][0],
															seq_feature.qualifiers["product"][0])
						else:
							prot_seq.description = "%s|%s" % (seq_record.name,
															seq_feature.qualifiers["product"][0])	
						entries.append(prot_seq)
	
	with open(options.output, "w") as handle:
		SeqIO.write(entries, handle, "fasta")
		
	
		
if __name__ == "__main__":
	__main__()
	