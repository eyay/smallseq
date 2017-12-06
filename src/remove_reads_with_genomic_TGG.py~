import dr_tools, sys, pysam, GenomeFetch
import os, argparse
sys.path.append('/mnt/crick/sandberglab/sra')
from joblib import Parallel, delayed

"""
This script removes reads when 3' has TGG on the genome
Works for smallrna star folder. aka max32, max38

>adapter_three_prime
TGGAATTCTCGGGTGCCAAGG
>polyA
AAAAAAAAAAAAA
"""

def safe_mkdir(path):
	if not os.path.exists(path):
		os.mkdir(path)
		os.chmod(path, 0o774)

def remove_reads_from_precursor(inbam):

	"""
	prepare input/output files
	"""
	inbamPysamObj = pysam.Samfile(inbam, "rb" )
	p = inbam.split("/")
	outbamTmp = "/".join(p[:-3]+[o.outstardir]+p[-2:])
	bam_out = ".".join(outbamTmp.split(".")[:-1]) + "_tmp.bam"
	bam_out_sorted = ".".join(outbamTmp.split(".")[:-1])
	outbam = pysam.Samfile(bam_out, "wb", template=inbamPysamObj)

	"""
	create genome fetch object
	"""
	gf = GenomeFetch.GenomeFetch(genomedir=o.genome_dir)

	"""
	remove reads when 3' has TGG on the genome
	"""
	for read in inbamPysamObj:
		read_name = read.qname
		tid = read.rname
		readchr  = inbamPysamObj.getrname(tid)
		readstart = int(read.pos) + 1
		readend = read.aend
		strand = read.flag
		readlen = len(read.seq) #this is the actual read length (41M, means readlen=41)
		read_len = read.qlen  #this only considers matches (8S30M, means read_len=30)
		minRlen = o.minRlen
		if readlen <= o.readlen_cutoff:
			outbam.write(read)
			continue
		
		if strand ==0:  #read maps to forward strand
			upperlimit = minRlen - readlen
			bpwindow = gf.get_seq_from_to(readchr, readend+1, readend+upperlimit)
			if readlen==minRlen-1 and (bpwindow == "T" or bpwindow == "A"): continue #TGGAATTCTCGGGTGCCAAGG
			elif readlen==minRlen-2 and (bpwindow == "TG" or bpwindow == "AA"): continue
			elif readlen==minRlen-3 and (bpwindow == "TGG" or bpwindow == "AAA"): continue
			elif readlen==minRlen-4 and (bpwindow == "TGGA" or bpwindow == "AAAA"): continue
			elif readlen==minRlen-5 and (bpwindow == "TGGAA" or bpwindow == "AAAAA"): continue
			else: outbam.write(read)

		elif strand ==16:  #read maps to reverse strand
			upperlimit = minRlen - readlen
			bpwindow = gf.get_seq_from_to(readchr, readstart-upperlimit, readstart-1)
			if readlen==minRlen-1 and (bpwindow == "A" or bpwindow == "T"): continue #TTCCA
			elif readlen==minRlen-2 and (bpwindow == "CA" or bpwindow == "TT"): continue
			elif readlen==minRlen-3 and (bpwindow == "CCA" or bpwindow == "TTT"): continue
			elif readlen==minRlen-4 and (bpwindow == "TCCA" or bpwindow == "TTTT"): continue
			elif readlen==minRlen-5 and (bpwindow == "TTCCA" or bpwindow == "TTTTT"): continue
			else: outbam.write(read)

	outbam.close()
	#sort and index the final bam file
	pysam.sort(bam_out, bam_out_sorted)	
	pysam.index(bam_out_sorted+".bam", template=inbamPysamObj)
	os.remove(bam_out)

#main function
if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--instardir', required=True)
	parser.add_argument('-o', '--outstardir', required=True)
	parser.add_argument('-g', '--genome_dir', default="/mnt/kauffman/ilqara/genomes/hg38_gtrnadb/")
	parser.add_argument('-c', '--readlen_cutoff', default=35)
	parser.add_argument('-x', '--minRlen', default=41, type=int)  #minimum read length to define a precursor
	parser.add_argument('-p', '--numCPU', default=20, type=int)
	o = parser.parse_args()

	if not os.path.exists(o.outstardir): safe_mkdir(o.outstardir)

	sample_names = os.listdir(o.instardir)
	samplenames_with_fullpath = []
	for sample in sample_names:
		##prepare input files
		bam = os.path.join(o.instardir, sample, "%s.bam" %sample)
		samplenames_with_fullpath.append(bam)

		path_outbam = os.path.join(o.outstardir, sample)
		if not os.path.exists(path_outbam): safe_mkdir(path_outbam)

	Parallel(n_jobs=o.numCPU)(delayed(remove_reads_from_precursor)(sample) for sample in samplenames_with_fullpath)

