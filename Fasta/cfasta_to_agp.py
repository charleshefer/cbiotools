###############################################################################
#Create an AGP file from a fasta file
#
#@author:charles.hefer@agresearch.co.nz
#@version:0.1
###############################################################################
import optparse
from Bio import SeqIO


def save_entry(object, object_beg, object_end, part_number, component_type = "N"):
	#reset the object_end
	
	print("%s\t%s\t%s\t%s\t%s" % (object, object_beg, object_end, part_number, component_type))
	
	object_beg = object_end + 1
	object_end = object_beg
	part_number = part_number + 1
	return(object_beg, object_end, part_number)


def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input file")
	parser.add_option("-o", "--ouput", default=None, dest="output",
					  help="The output .agp file")
	(options, args) = parser.parse_args()
	
	if not options.input:
		parser.error("Need to specify the input file")
	if not options.output:
		parser.error("Need to specify the output directory")
		
	nucleotides = ["A","C","T","G"]
	space = ["N"]
	nucleotide_run = None
	space_run = None
	
	with open(options.input,"r") as handle:
		for entry in SeqIO.parse(handle, "fasta"):
			object_beg = 1
			object_end = 1
			part_number = 1
			nucleotide_run = None
			print(str(entry.seq))
			for base in str(entry.seq):
				if base in nucleotides and nucleotide_run is None:
					#this is the first nucleotide sequence
					object_end = object_end + 1				
					nucleotide_run = True
					continue

				if base in space and nucleotide_run is None:
					#this is the first space run
					object_end = object_end + 1
					nucleotide_run = False
					continue

				
				if base in nucleotides and nucleotide_run is True:
					#continue the nucleotide run
					object_end = object_end + 1
					continue

				if base in space and nucleotide_run is False:
					#continue the space run
					object_end = object_end + 1
					continue

				if base in space and nucleotide_run is True:
					#this is a change from a nucleotide run to a space run
					#the current base is the end of a nucleotide run, subtract one
					object_end = object_end - 1
					# save the entry (previous section)
					# the previous section was a wgs component
					component_type = "W"					
					object_beg, object_end, part_number = save_entry(entry.id, object_beg, object_end, part_number, component_type)
					#reset to the current position
					#object_beg = object_end
					object_end = object_end + 1
					
					#is not a nucleotide run any more
					nucleotide_run = False
					continue

				if base in nucleotides and nucleotide_run is False:
					#Change from a space run to a nubleotide run
					#reset the end position with one
					object_end = object_end - 1
					#save the entry
					object_beg, object_end, part_number = save_entry(entry.id, object_beg, object_end, part_number)
					#reset to the current position
					#object_beg = object_end
					object_end = object_end + 1
					
					nucleotide_run = True
					continue
			save_entry(entry.id, object_beg, object_end, part_number)
			

if __name__ == "__main__":
	__main__()
	