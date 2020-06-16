###############################################################################
#Download a file with sequence entries from the NCBI using EUTILS
#@requires:biopython
#@author:charles.hefer@agresearch.co.nz
#@version:0.1
###############################################################################
import optparse
from Bio import Entrez, SeqIO
import time, os, shutil
import tempfile


Entrez.email="charles.hefer@agresearch.co.nz"
API_KEY = "fbe3971629fd3304f0ce1c599121182d1007"

def __main__():
    """Parse the cmd lne options"""
    parser = optparse.OptionParser()
    parser.add_option("-f", "--file", default=None, dest="file",
    				  help="The file containint accessions accession to retrieve")
    parser.add_option("-d", "--database", default=None, dest="database",
					help = "The database to retrieve from")
    parser.add_option("-o", "--output", default=None, dest="output",
					  help="The output file")
    (options, args) = parser.parse_args()
	
    if not options.file:
        parser.error("Need to specify the input file")
    if not options.database:
        parser.error("Need to specify the database")
    if not options.output:
        parser.error("Need to specify the output file")

    seqrecords = []

    with open(options.file) as fhandle:
        for line in fhandle:
            line = line.rstrip()

            if line.startswith("gi"):
                try:
                    handle = Entrez.efetch(db=options.database,
                        id = line,
                        rettype="fasta",
                        api_key=API_KEY)
                    time.sleep(0.5)
                    
                #This is not a fasta object.... it is just text
                #try:
                    
                    record = handle.read()
                    
                
                    #Save the record to file
                    fasta_handle = tempfile.NamedTemporaryFile(delete=False, mode="w")
                    fasta_handle.write(record)
                    fasta_handle.close()

                    #  Now, modify the record, add GI to the identifiyer
                    with open(fasta_handle.name, "r") as ncbi_handle:
                        for record in SeqIO.parse(ncbi_handle, "fasta"):
                            #print(record.description)
                            record.description = record.description.replace(record.id, "")
                            record.id = record.id + "|" + line
                        
                            seqrecords.append(record)
                except Exception as e:
                    print("***********%s NOT FOUND**************" % line)
                    
    with open(options.output, "w") as outhandle:
        SeqIO.write(seqrecords, outhandle, "fasta")

if __name__ == "__main__":
    __main__()
	