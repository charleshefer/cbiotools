###############################################################################
#cbam_deepest_coverage.py
#Given a bam file (and index) and a list of chromosomes, this returns the 
# position with the highest coverage
# Uses samtools depth
#
#@author:charles.hefer@agresearch.co.nz
#@version:0.1
###############################################################################
import optparse
import subprocess


def __main__():
	"""Parse the cmd lne options"""
	parser = optparse.OptionParser()
	parser.add_option("-i", "--input", default=None, dest="input",
					  help="The input bam file")
	parser.add_option("-c", "--chromfile", default=None, dest="chromfile",
					  help="The input list of chromosomes ")
	parser.add_option("-d", "--directory", default=None, dest="directory",
						help="directory")
	parser.add_option("-l", "--lane", default=None, dest="lane",
						help="Lane")				  
	parser.add_option("-o", "--output", default=None, dest="output",
					  help="The output file")
	(options, args) = parser.parse_args()
	
	if not options.input:
		parser.error("Need to specify the input file")
	if not options.chromfile:
		parser.error("Need to specify the list of chromosomes file")
	if not options.output:
		parser.error("Need to specify the output file")

	with open(options.output, "w") as outfile:
		with open(options.chromfile, "r") as handle:
			for line in handle:
				line = line.rstrip()
				#p = subprocess.run(["echo %s" % line], shell=True, check=True, stdout=subprocess.PIPE)
				p = subprocess.run(["samtools depth -r %s %s/%s | sort -n -k 3 | tail -n 1" % (line, options.directory, options.input)], shell=True, check=True, stdout=subprocess.PIPE)
				outfile.write(p.stdout.decode('utf-8').rstrip() + "\t%s" % options.input + "\t%s\n" % options.lane)

		
if __name__ == "__main__":
	__main__()
	