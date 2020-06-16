###############################################################################
#Download a file with sequence entries from the NCBI using EUTILS
#@requires:biopython
#@author:charles.hefer@agresearch.co.nz
#@version:0.1
###############################################################################
import optparse
from Bio import Entrez, SeqIO
import time, os


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

    try:
        os.mkdir("./tmp/")
    except OSError as e:
        pass


    with open(options.file) as fhandle:
        for line in fhandle:
            line = line.rstrip()

            if line.startswith("gi"):
                handle = Entrez.efetch(db=options.database,
                    id = line,
                    rettype="fasta",
                    api_key=API_KEY)
                time.sleep(0.3)

                record = handle.read()
                #This is not a fasta object.... it is just text
            
                #Save the record to file
                with open("./tmp/" + line + ".fasta", "w") as outhandle:
                    outhandle.write(record)

                #  Now, modify the record, add GI to the identifiyer
                with open("./tmp/" + line + ".fasta", "r") as ncbi_handle:
                    for record in SeqIO.parse(ncbi_handle, "fasta"):
                        #there should be only one record
                        record.id = record.id + "|" + line
                        #record.description = record.description
                        with open("./tmp/" + line + "_mod.fasta", "w") as record_handle:
                            SeqIO.write(record, record_handle, "fasta")
                    

    seqrecords = []

    #read all the records from the tempfile
    for file in os.listdir("./tmp/"):
        if "_mod" in file: 
            handle = open("./tmp/" + file, "r")
            for record in SeqIO.parse(handle, "fasta"):
                seqrecords.append(record)

    with open(options.output, "w") as outhandle:
        SeqIO.write(seqrecords, outhandle, "fasta")

if __name__ == "__main__":
    __main__()
	