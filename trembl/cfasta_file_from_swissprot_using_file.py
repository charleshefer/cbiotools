###############################################################################
#Download a file with sequence entries from the Swissprot using ExPASY
#@requires:biopython
#@author:charles.hefer@agresearch.co.nz
#@version:0.1
###############################################################################
import optparse
#Not using the SwissProt parser, issues with the new format released 2020-05-01
from Bio import SeqIO
import time, os, shutil
import urllib
import tempfile

def __main__():
    """Parse the cmd lne options"""
    parser = optparse.OptionParser()
    parser.add_option("-f", "--file", default=None, dest="file",
    				  help="The file containing accession to retrieve")
    parser.add_option("-o", "--output", default=None, dest="output",
					  help="The output file")
    (options, args) = parser.parse_args()
	
    if not options.file:
        parser.error("Need to specify the input file")
    if not options.output:
        parser.error("Need to specify the output file")


    seqrecords = []
    with open(options.file) as fhandle:
        for line in fhandle:
            line = line.rstrip()

            #This usecase contained "tr|" followed by accession as the line

            #Only look at trembl entries
            if line.startswith("tr|"):
                accession = line.split("|")[-1]
                
                #grab the xml entry from uniprot
                with urllib.request.urlopen("https://www.uniprot.org/uniprot/%s.xml" % accession) as response: 
                    record = response.read()

                #Save the record to a temp file
                xml_record = tempfile.NamedTemporaryFile(delete=False)
                xml_record.write(record)
                xml_record.close()

                #Convert the xml file to a fasta file
                with open(xml_record.name, "r") as xml_handle:
                    for record in SeqIO.parse(xml_handle, "uniprot-xml"):
                        record.id = record.id + "|" + line
                        seqrecords.append(record)

    with open(options.output, "w") as outhandle:
        SeqIO.write(seqrecords, outhandle, "fasta")

if __name__ == "__main__":
    __main__()
	