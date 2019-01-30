###############################################################################
#Convert genbank to fasta format
#
#@requires Biopython
#@author:charles.hefer@gmail.com
#@version:0.1b
###############################################################################
import optparse, sys
from Bio import SeqIO
import copy

def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input file")
	parser.add_option("-p", "--ptt", default=None, dest="ptt",
					  help="The output ptt file")
	parser.add_option("-r", "--rnt", default=None, dest="rnt",
					  help="The output rnt file")
	
	(options, args) = parser.parse_args()
	
	if not options.input:
		parser.error("Need to specify the input genbank file")
	if not options.ptt:
		parser.error("Need to specify the output ptt file")
	if not options.rnt:
		parser.error("Need to specify the output rnt file")
		
	proteins = []
	rnas = []
	
	header = {}
	
	coldetail={"Location":"",
			  "Strand":"",
			  "Length":"",
			  "PID":"",
			  "Gene":"",
			  "Synonym":"",
			  "Code":"",
			  "COG":"",
			  "Product":""}
	
	with open(options.input, "r") as handle:
		for entry in SeqIO.parse(handle, "genbank"):
			#the header
			header["description"] = entry.description
			header["proteins"] = 0
			header["rna"] = 0
			for feature in entry.features:
				if "translation" in feature.qualifiers.keys():
					#print(feature.__dict__)
					header["proteins"] = header["proteins"] + 1
					
					#Figure out the strandedness and save the location(position)
					if feature.location.strand > 0:
						coldetail["Location"] = "%i..%i" % (feature.location._start,
													   feature.location._end)
						coldetail["Strand"] = "+"
					else:
						coldetail["Location"] = "%i..%i" % (feature.location._end,
													   feature.location._start)
						coldetail["Strand"] = "-"
						
					#Calculate the feature length
					coldetail["Length"] = str(int(feature.location._end) -
												  int(feature.location._start))
					
					#PID
					try:
						coldetail["PID"] = feature.qualifiers["protein_id"][0]
					except KeyError:
						coldetail["PID"] = feature.qualifiers["locus_tag"][0]
					#Gene
					try:
						coldetail["Gene"] = feature.qualifiers["gene"][0]
					except KeyError:
						coldetail["Gene"] = "-"
					
					#Synonym
					coldetail["Synonym"] = feature.qualifiers["locus_tag"][0]
					
					#Code
					coldetail["Code"] = "-"
					
					#COG
					coldetail["COG"] = "-"
					
					#Product
					coldetail["Product"] = feature.qualifiers["product"][0]
					
					proteins.append(copy.deepcopy(coldetail))
					#print("\n")
				
				if feature.type == "tRNA" or feature.type == "rRNA":
					#working with RNA here
					header["rna"] = header["rna"] + 1
					
					#Figure out the strandedness and save the location(position)
					if feature.location.strand > 0:
						coldetail["Location"] = "%i..%i" % (feature.location._start,
													   feature.location._end)
						coldetail["Strand"] = "+"
					else:
						coldetail["Location"] = "%i..%i" % (feature.location._end,
													   feature.location._start)
						coldetail["Strand"] = "-"
						
					#Calculate the feature length
					coldetail["Length"] = str(int(feature.location._end) -
												  int(feature.location._start))
					
					#PID
					try:
						coldetail["PID"] = feature.qualifiers["protein_id"][0]
					except KeyError:
						coldetail["PID"] = feature.qualifiers["locus_tag"][0]
					#Gene
					try:
						coldetail["Gene"] = feature.qualifiers["gene"][0]
					except KeyError:
						coldetail["Gene"] = "-"
					
					#Synonym
					coldetail["Synonym"] = feature.qualifiers["locus_tag"][0]
					
					#Code
					coldetail["Code"] = "-"
					
					#COG
					coldetail["COG"] = "-"
					
					#Product
					coldetail["Product"] = feature.qualifiers["product"][0]
					
					#print("\n")
					rnas.append(copy.deepcopy(coldetail))
					
				
	
	with open(options.ptt, "w") as handle:
		handle.write(header["description"] + "\n")
		handle.write("%s proteins\n" % header["proteins"])
		#header_line = [key for key in entry.keys()]
		handle.write("Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n")
		for entry in proteins:
			handle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (entry["Location"],
																   entry["Strand"],
																   entry["Length"],
																   entry["PID"],
																   entry["Gene"],
																   entry["Synonym"],
																   entry["Code"],
																   entry["COG"],
																   entry["Product"]))
					
			
	with open(options.rnt, "w") as handle:
		handle.write(header["description"] + "\n")
		handle.write("%s RNAs\n" % header["rna"])
		#header_line = [key for key in entry.keys()]
		handle.write("Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n")
		for entry in rnas:
			handle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (entry["Location"],
																   entry["Strand"],
																   entry["Length"],
																   entry["PID"],
																   entry["Gene"],
																   entry["Synonym"],
																   entry["Code"],
																   entry["COG"],
																   entry["Product"]))
		
		
	
		
if __name__ == "__main__":
	__main__()
	