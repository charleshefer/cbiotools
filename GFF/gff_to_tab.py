###############################################################################
#Converts gff to tab
#
#@author:charles.hefer@agresearch.co.nz
#@version:0.1
###############################################################################
import optparse
import re


id_re = re.compile("ID=(?P<ID>.+?);")
name_re = re.compile(".+?Name=(?P<name>.+?);")
locus_tag_re = re.compile(".+?locus_tag=(?P<locus_tag>.+?);")
old_locus_tag_re = re.compile(".+?old_locus_tag=(?P<old_locus_tag>.+?);")
protein_id_re = re.compile(".+?protein_id=(?P<protein_id>.+?);")
product_re = re.compile(".+?product=(?P<product>.+?);")


def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input gff file")
	parser.add_option("-o", "--ouput", default=None, dest="output",
					  help="The output tab file")
	(options, args) = parser.parse_args()
	
	if not options.input:
		parser.error("Need to specify the input file")
	if not options.output:
		parser.error("Need to specify the output file")


	outlines = []

	with open(options.input, "r") as handle:
		for line in handle:

			id = "NA"
			id_description = "NA"
			name = "NA"
			locus_tag = "NA"
			old_locus_tag = "NA"
			protein_id = "NA"
			product = "NA"


			if line.startswith("#"):
				continue
			#print(line)
			cols = line.split("\t")
			id = cols[0]
			description = cols[8].rstrip() + ";"
			
			if id_re.match(description):
				id_description = id_re.match(description)["ID"]
			if name_re.match(description):
				name = name_re.match(description)["name"]
			if locus_tag_re.match(description):
				locus_tag = locus_tag_re.match(description)["locus_tag"]
			if old_locus_tag_re.match(description):
				old_locus_tag = old_locus_tag_re.match(description)["old_locus_tag"]
			if protein_id_re.match(description):
				protein_id = protein_id_re.match(description)["protein_id"]
			if product_re.match(description):
				product = product_re.match(description)["product"]
			
			outlines.append("\t".join([id, id_description,name, 
			locus_tag, old_locus_tag, protein_id,
			product + "\n"]))

	with open(options.output, "w") as outhandle:
		outhandle.write("\t".join(["id", "id_description", "name", "locus_tag", "old_locus_tag", 
						"protein_id","product" + "\n"]))
		for outline in outlines:
			outhandle.write(outline)
		
if __name__ == "__main__":
	__main__()
	