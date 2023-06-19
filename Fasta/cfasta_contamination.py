###############################################################################
#Given a contamination file from genbank, parses and fix the issues
#
#@requires: biopython
#@author:charles.hefer@gmail.com
#@version:0.1
###############################################################################
import optparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import copy

def __main__():
	"""Parse the cmd line options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input file")
	parser.add_option("-o", "--ouput", default=None, dest="output",
					  help="The output file")
	parser.add_option("-c", "--contamination", default=None, dest="contamination",
					help = "The contamination file")
	(options, args) = parser.parse_args()
	
	
	if not options.input:
		parser.error("Need to specify the input file")
	if not options.output:
		parser.error("Need to specify the output file")
	if not options.contamination:
		parser.error("Need to specify a contamination file")

	exclude_entries = []
	mask_entries = []
	mask_entry_list = []
	mask_entry = {"id" : None,
								"positions" : []}

	with open(options.contamination, "r") as handle:
		exclude = False
		exclude_list = False
		trim = False
		for line in handle:
			line = line.rstrip()
			if line.startswith("Exclude:"): 
				exclude = True
				continue
			if exclude:
				cols = line.split("\t")
				#these are the entries that needs to be excluded
				if len(cols) == 3:
					exclude_list = True
					exclude_entries.append(cols[0])
					
			if exclude and exclude_list:
					cols = line.split("\t")
					if len(cols) != 3:
						exclude = False
						exclude_list = False
			
			if line.startswith("Trim:"):
				trim = True
				continue
			
			if trim:
				cols = line.split("\t")
				
				if len(cols) == 4:
					positions = cols[2]
					mask_entry["id"] = cols[0]
					mask_entry["positions"] = []
					for position in positions.split(","):
						start = int(position.split("..")[0]) - 1
						stop = int(position.split("..")[1]) - 1
						#mask_entry["positions"]["start"] = int(position.split("..")[0]) - 1
						#mask_entry["positions"]["end"] = int(position.split("..")[1]) - 1
						mask_entry["positions"].append([start, stop])
					mask_entries.append(copy.deepcopy(mask_entry))
					mask_entry_list.append(cols[0])
	
	fasta_entries = []
	with open(options.input, "r") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			entry = {"id" : record.id,
					"description" : "",
					"seq" : str(record.seq)}
			fasta_entries.append(entry)

	for entry in fasta_entries:
		if entry["id"] not in exclude_entries:
			if entry["id"] in mask_entry_list:
				for mask_entry in mask_entries:
					if mask_entry["id"] == entry["id"]:
						for position in mask_entry["positions"]:
							print("Before replacement:")
							print(len(entry["seq"]))
							length_of_ns = "N" * (position[1] - position[0])
							entry["seq"] = entry["seq"][:position[0]] + length_of_ns + entry["seq"][position[1]:]
							print("After replacement:")
							print(len(entry["seq"]))

	fasta_records = []
	for entry in fasta_entries:
		new_record = SeqRecord(Seq(entry["seq"]))
		new_record.id = entry["id"]
		new_record.description = entry["description"]
		fasta_records.append(new_record)

	
	with open(options.output, "w") as outhandle:
		SeqIO.write(fasta_records, outhandle, "fasta")
	
if __name__ == "__main__":
	__main__()