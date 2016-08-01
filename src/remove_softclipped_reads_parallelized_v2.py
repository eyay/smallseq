from __future__ import division
import dr_tools, sys, os, pysam, argparse
from joblib import Parallel, delayed

"""
This script:
	- removes all hardclipped reads, reads with insertions and deletions
	- removes reads with any 5' soft-clipping
	- removes reads with more than 3 bases on 3' via soft-clipping 
 
"""

def safe_mkdir(path):
	if not os.path.exists(path):
		os.mkdir(path)
		os.chmod(path, 0o774)

def remove_clipped_reads(inbam):
	inbamPysamObj = pysam.Samfile(inbam, "rb" )
	p = inbam.split("/")
	outbamTmp = "/".join(p[:-3]+[o.outstardir]+p[-2:])
	bam_out = ".".join(outbamTmp.split(".")[:-1]) + "_tmp.bam"
	bam_out_sorted = ".".join(outbamTmp.split(".")[:-1])
	outbam = pysam.Samfile(bam_out, "wb", template=inbamPysamObj)
	total_reads = 0; reads_after_clip_removal = 0
	for read in inbamPysamObj:
		total_reads += 1
		read_name = read.qname
		tid = read.rname
		readchr  = inbamPysamObj.getrname(tid)
		cigar = read.cigar  #this is a list of tuples;  example: cigar=[(0, 25), (4, 8)]  

		num_softclipped_5p = 0
		num_softclipped_3p = 0
		num_hardclipped = 0
		num_insertions = 0
		num_deletions = 0

		for i,e in enumerate(cigar):
			if e[0] == 5:
				num_hardclipped += e[1]
			elif e[0] == 1:
				num_insertions += e[1]
			elif e[0] == 2:
				num_deletions += e[1]
			elif i == 0 and e[0] == 4:
				num_softclipped_5p += e[1]
			elif i == len(cigar)-1 and e[0] == 4:
				num_softclipped_3p += e[1]

		if num_softclipped_5p > o.allowed_num_clipped_bases_5p: continue
		elif num_softclipped_3p > o.allowed_num_clipped_bases_3p: continue
		elif num_hardclipped > 0: continue
		elif num_insertions > 0: continue 
		elif num_deletions > 0: continue   

		outbam.write(read)
		reads_after_clip_removal += 1
	outbam.close()
	#sort and index the final bam file
	pysam.sort(bam_out, bam_out_sorted)     
	pysam.index(bam_out_sorted+".bam", template=inbamPysamObj)

	print reads_after_clip_removal, total_reads
	print "%s percentage of removed reads = %.2f" %(p[-1], 100*(1 -reads_after_clip_removal/total_reads))

	os.remove(bam_out)

#main function
if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--instardir', required=True)
	parser.add_argument('-o', '--outstardir', required=True)
	parser.add_argument('-n', '--allowed_num_clipped_bases_5p', type=int, default=0, help='both hard and soft clipped, S:4, H:5')
	parser.add_argument('-t', '--allowed_num_clipped_bases_3p', type=int, default=3, help='both hard and soft clipped, S:4, H:5')
	parser.add_argument('-p', '--numCPU', default=20, type=int)
	o = parser.parse_args()

	if not os.path.exists(o.outstardir): safe_mkdir(o.outstardir)

	sample_names = os.listdir(o.instardir)
	samplenames_with_fullpath = []
	for sample in sample_names:
		##prepare input files
		unique_bam = os.path.join(o.instardir, sample, "%s_unique.bam" %sample)
		multi_bam = os.path.join(o.instardir, sample, "%s_multi.bam" %sample)
		samplenames_with_fullpath.append(unique_bam)
		samplenames_with_fullpath.append(multi_bam)

		path_outbam = os.path.join(o.outstardir, sample)
		if not os.path.exists(path_outbam): safe_mkdir(path_outbam)

	Parallel(n_jobs=o.numCPU)(delayed(remove_clipped_reads)(sample) for sample in samplenames_with_fullpath)

