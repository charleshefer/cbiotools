###############################################################################
#Fix the Ns
#For Andrew Griffiths
#@author:charles.hefer@agresearch.co.nz
#@version:0.1
###############################################################################
import optparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

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

	space_re = re.compile(r"[A,C,G,T]*(N+)[A,C,G,T]*")
	replace_re = re.compile(r"[A,C,G,T]+?(N{2})[A,C,G,T]+?")
	records = []


	with open(options.input, "r") as handle:
		for entry in SeqIO.parse(handle, "fasta"):
			if entry.id == "chr0":
				#split all the NNs from the sequence, no matter how many
				nlength = []
				res = space_re.findall(str(entry.seq))
				if res:
					for r in res:
						nlength.append(len(r))
					#make unique and sort in reverse order
					nlength = list(set(nlength))
					nlength.reverse()
					seq = str(entry.seq)
					
					#replace from largest to smallest
					for nlen in nlength:
						seq = seq.replace("N"*nlen, "--CUT-HERE--")
					
					#split, rename and append to the records
					counter = 1
					for s in seq.split("--CUT-HERE--"):
						if len(s) > 0:
							sequence = SeqRecord(Seq(s))
							sequence.id = "chr0_" + str(counter)
							sequence.description = ""
							counter = counter + 1 
							records.append(sequence)

			else:
				#only replace the occurance of exactly 10000 Ns
				#print(entry.id)
				#print(entry.seq)
				res = replace_re.search(str(entry.seq))
				seq = str(entry.seq)
				while res:
					seq = seq[:res.span(1)[0]] + "N"*100 + seq[res.span(1)[1]:]
					res = replace_re.search(seq)
				#print(seq)
				sequence = SeqRecord(Seq(seq))
				sequence.id = entry.id
				sequence.description = ""
				records.append(sequence)
				#for i in range(res.span(1)[0],res.span(1)[1]):
				#	print(str(entry.seq)[i])
				


	with open(options.output, "w") as handle:
		SeqIO.write(records, handle, "fasta")




if __name__ == "__main__":
	__main__()
	