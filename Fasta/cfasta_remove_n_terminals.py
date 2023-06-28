###############################################################################
#cfasta_remove_n_terminals
#Reads in a fasta formatted file, replaces NNs in the beginning of a sequence
# with AAs, and removes any NNs from the end of the entries
#
#@requires biopython
#@author:charles.hefer@gmail.com
#@version:0.1
###############################################################################
import optparse
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random


def random_base(n):
		"""Returns a string of length n"""
		seq = []
		for i in range(0,n):
			bases = ["A","C","G","T"]
			n = random.randint(0,3)
			seq.append(bases[n])
		return "".join(seq)



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
	
	start_with_N = re.compile("^(N+?)[A,C,T,G]")
	end_with_N = re.compile("^.+?(N+?)$")
	
	fasta_records = []
	start_with_NNN = False
	end_with_NNN = False

	

	
	with open(options.input, "r") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			print(record.id)
			seq = str(record.seq)
			start_with_NNN = False
			end_with_NNN = False

			if start_with_N.match(seq):
				start_with_NNN = True
				start_groups = start_with_N.match(seq).groups()[0]
			
			if end_with_N.match(seq):
				end_with_NNN = True
				end_groups = end_with_N.match(seq).groups()[0]
			
			if start_with_NNN or end_with_NNN:

				if start_with_NNN:
					#print(len(seq))
					seq = random_base(len(start_groups)) + seq[len(start_groups):]
					#print(len(seq))
				
				if end_with_NNN:
					seq = seq[:-len(end_groups)]
			
				if len(seq) > 1:
					sequence = SeqRecord(Seq(seq))
					sequence.id = record.id
					sequence.description = ""
					fasta_records.append(sequence)
				else:
					break
			else:
				fasta_records.append(record)
	
	with open(options.output, "w") as outhandle:
		SeqIO.write(fasta_records, outhandle, "fasta")
		
if __name__ == "__main__":
	__main__()
	