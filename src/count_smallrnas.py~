from __future__ import division
from collections import defaultdict
import dr_tools, pysam
import argparse, os
from joblib import Parallel, delayed


def safe_mkdir(path):
	if not os.path.exists(path):
		os.mkdir(path)
		os.chmod(path, 0o774)

class betweenRE(dr_tools.Cregion):
	pass


def bam_to_windows(inbam):
	inbamPysamObj = pysam.Samfile(inbam, "rb" )
	p = inbam.split("/")
	samplename = p[-2]
	outbamTmp = "/".join(p[:-3]+[o.outdir]+p[-2:])
	tempCountfile = ".".join(outbamTmp.split(".")[:-1]) + "_tmpCount.txt"
	finalCountfile = ".".join(outbamTmp.split(".")[:-1]) + "_Count.txt"
	read2overlapCoords=defaultdict(list)

	for read in inbamPysamObj:
		readname = read.qname
		tid = read.rname
		readchr  = inbamPysamObj.getrname(tid)
		readstart = int(read.pos)
		readend = read.aend
		if read.is_reverse: 
			strand="-"
		else:
			strand="+"
		readlen = len(read.seq) #this is the actual read length (41M, means readlen=41)
		read_len = read.qlen  #this only considers matches (8S30M, means read_len=30)
 
		midpos = (readstart + readend)//2

		#retrieve list of overlapping coordinates
		overlap_list = betweenRE.overlappingpoint(readchr, midpos, strand)
		annotatedCount = len(overlap_list)

		#make a dictionary of read and overlapping coordinates
		read2overlapCoords[readname].append(overlap_list)

	with open(tempCountfile, "w") as outfh:
		for read in read2overlapCoords:
			coordsList = read2overlapCoords[read]
			readCount = len(coordsList)
			annotatedCount = readCount-coordsList.count([])
			#len(coordsList) is never zero
			for coord in coordsList:
				if len(coord) == 0:
					print >> outfh, dr_tools.join(read, "NA", readCount, annotatedCount)
				else:
					###coord[1] will be double-counting
					coord = str(coord[0])  #otherwise I got keyError. it was "instance" type variable
					geneid = coord2geneid.get(coord, 'NA')
					print >> outfh, dr_tools.join(read, geneid, readCount, annotatedCount)
	outfh.close()

	## readCount, annotatedCount scenarios
	# 1, 1  unique map, annotated to single gene, counts as 1
	# 2, 1  multi map, annotated to single gene, count as 1, discard other alignment
	# n, n  multi map, annotated to two genes, count as 1/n
	# k, m  where k>m and m>1, multi map, annotated to multi genes, count 1/m, discard other alignment
	#
	#formula is always: count = 1/annotatedCount
	geneid2counts={}
	unannotReadsDict={}
	for line in open(tempCountfile, "r"):
		p = line.split()
		read, geneid, readCount, annotatedCount = p
		annotatedCount = int(annotatedCount)
		if not geneid in geneid2counts: geneid2counts[geneid] = 0
		if annotatedCount > 0:
			geneid2counts[geneid] += 1/annotatedCount
		else:
			geneid2counts[geneid] += 0
		if annotatedCount < readCount and annotatedCount == 0:
			unannotReadsDict[read] = 1

	num_unannot = len(unannotReadsDict)
	num_annot = 0
	for geneid in geneid2counts:
		if "P-cel" in geneid: continue
		if geneid == "NA": continue
		num_annot += geneid2counts[geneid]

	with open(finalCountfile, "w") as outfh2:
		print >> outfh2, dr_tools.join("#samples", samplename)
		print >> outfh2, dr_tools.join("#unannotatedmolc", num_unannot)
		print >> outfh2, dr_tools.join("#annotatedmolc", num_annot)
		for geneid in geneidlist:
			print >> outfh2, dr_tools.join(geneid2name[geneid], geneid, geneid2counts.get(geneid, "0"))
	outfh2.close()

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--instardir', required=True)
	parser.add_argument('-o', '--outdir', required=True, help='directory for final counts')
	parser.add_argument('-g', '--GenePred', default="annotations/combined_annots.gp") 
	parser.add_argument('-p', '--numCPU', default=20, type=int)
	o = parser.parse_args()

	### read coordinates and add to windows
	geneid2name={}
	coord2geneid={}
	geneidlist=[]
	for line in open(o.GenePred, 'r'):
		p = line.split()
		chrom, startpos, endpos, strand, geneid, genename = p[2], int(p[4]), int(p[5]), p[3], p[1], p[12]
		#"add to windows
		r = betweenRE(chrom, startpos, endpos, strand)
		r.addtowindows()
		coord = chrom+":"+str(startpos+1)+"-"+str(endpos)+":"+strand
		coord2geneid[coord] = geneid #check length
		geneid2name[geneid] = genename
		geneidlist.append(geneid)

	### create outfolder, add sample names to list
	if not os.path.exists(o.outdir): safe_mkdir(o.outdir)
	sample_names = os.listdir(o.instardir)
	samplenames_with_fullpath = []
	for sample in sample_names:
		##prepare input files
		uniqmultibam = os.path.join(o.instardir, sample, "%s.bam" %sample)
		samplenames_with_fullpath.append(uniqmultibam)
		path_outbam = os.path.join(o.outdir, sample)
		if not os.path.exists(path_outbam): safe_mkdir(path_outbam)

	### call function in parallel 
	Parallel(n_jobs=o.numCPU)(delayed(bam_to_windows)(sample) for sample in samplenames_with_fullpath)

