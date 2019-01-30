# Working with genbank files
Genbank files (.gb, .gbk) can be large and quite compled. Here are some tools to extract information from these files into some common, and not so common formats.

## cgenbank_to_fasta.py
Converts a genbank file to a fasta file. This requires BioPython to be installed.

## cgenbank_to_ptt.py
Converts a genbank file to an ncbi protein (.ptt) and ncbi rna (.rnt) file. These are not widely used (or might be, but I have never seen this. Requires BioPython. This produces files that can be imported into `Rockhopper` (http://cs.wellesley.edu/~btjaden/Rockhopper/user_guide.html)

