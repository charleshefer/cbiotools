###############################################################################
#cfasta_rearanges
#Reads in a fasta formatted file, identifies a contig, cut a piece of sequence,
#and moves it to the end
#can also reverse the sequence
#
#@requires biopython
#@author:charles.hefer@gmail.com
#@version:0.1
###############################################################################
import optparse, sys
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
	parser.add_option("-s", "--start", default=None, dest="start",
					  help="The start position to cut")
	parser.add_option("-e", "--end", default=None, dest="end",
					  help="The end position to cut")
	parser.add_option("-p", "--paste", default=None, dest="paste",
					  help="The insert/paste position")
	parser.add_option("-c", "--contig", default=None, dest="contig",
					  help="The name of the contig to cut and paste from")
	parser.add_option("-r", "--reverse", default=False, dest="reverse",
					  help="Reverse the cut sequence before pasting")
	(options, args) = parser.parse_args()
	
	if not options.input:
		parser.error("Need to specify the input file")
	if not options.output:
		parser.error("Need to specify the output file")
	if not options.contig:
		parser.error("Need to specify the name of the contig to work with")
	if not options.start:
		parser.error("Need to specify cut start position")
	if options.start:
		options.start = int(options.start)
	if not options.end:
		parser.error("Need to specify cut end position")
	if options.end:
		options.end = int(options.end)
	if not options.paste:
		parser.error("Need to specify insert/paste position")
	if options.paste:
		options.paste = int(options.paste)

	#Read into RAM, will bomb if there are Gbs of entries
	#but this has not happened yet.
	#Better than IO overhead for most cases
	fasta_records = []
	
	with open(options.input, "r") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			if record.id == options.contig:
				print("Looking for contig to edit...")
				print("Found: %s" % record.id)
				#print(str(record._seq))
				record_len = len(str(record._seq))
				print("Record is %s bases long." % record_len)
				if options.start > record_len:
					print("Start cut position larger than record length: %s" % options.start)
					sys.exit()
				if options.end > record_len:
					print("Stop cut position larger than record length: %s" % options.stop)
					sys.exit()
				if options.paste > record_len:
					print("Cut sequence will be appended to the sequence")
					options.paste = record_len
				if options.paste == -1:
					options.paste = record_len
				print("Start: %s\t Stop: %s\t Paste: %s" % (options.start, options.end, options.paste))
				
				
				
				sequence = record._seq
				cut_sequence = sequence[options.start:options.end]
				print("Cut sequence: %s" % str(cut_sequence))
				
				inserted_sequence = "".join([str(sequence[0:options.paste]),
												str(cut_sequence),
											str(sequence[options.paste + 1:-1])
												])
				#print("Inserted_sequence: %s" % str(inserted_sequence))
				
				
				final_sequence = "".join([str(inserted_sequence[0:options.start]), str(inserted_sequence[options.end:])])
				print("Final sequence: %s %s" % (str(final_sequence), len(str(final_sequence))))

				new_record = SeqRecord(Seq(final_sequence), id=record.id, description="")



			fasta_records.append(new_record)
	
	with open(options.output, "w") as outhandle:
		SeqIO.write(fasta_records, outhandle, "fasta")
		
if __name__ == "__main__":
	__main__()
	