from __future__ import division
import os
import dr_tools, sys, argparse
import subprocess
#from joblib import Parallel, delayed

def safe_mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)
        os.chmod(path, 0o774)

def num_unique_reads(bam_f):
	cmd = 'samtools view -c -bq 255 %s' %bam_f
	nreads = subprocess.check_output(cmd, shell=True)[:-1]
	return nreads

def subsample_1M_unique(inbam, outbam, frac):
    cmd = 'samtools view -s %f -bq 255 -o %s %s' %(frac, outbam, inbam)
    subprocess.call(cmd, shell = True)


def subsample_by_readlength(bam_f, readlength_cutoff, out_file):
    if o.min_or_max_readlen == "max":  cmd = 'samtools view -h %s | awk \'length($10) <= %i || $1~"@"\' | samtools view -bS - > %s' %(bam_f, readlength_cutoff, out_file)
    if o.min_or_max_readlen == "min":  cmd = 'samtools view -h %s | awk \'length($10) >= %i || $1~"@"\' | samtools view -bS - > %s' %(bam_f, readlength_cutoff, out_file)
    subprocess.call(cmd, shell = True)


if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--stardir', required=True)
	parser.add_argument('-o', '--outdir', required=True)
	parser.add_argument('-c', '--readlen_cutoff', default=28, type=int)
	parser.add_argument('-m', '--min_or_max_readlen', required=True)
	o = parser.parse_args()

	if not os.path.exists(o.outdir): safe_mkdir(o.outdir)

	sample_names = os.listdir(o.stardir)
	for sample in sample_names:
		##prepare files and proper paths
		unique_bam = os.path.join(o.stardir, sample, "%s_unique.bam" %sample)
		multi_bam = os.path.join(o.stardir, sample, "%s_multi.bam" %sample)
		path_newbam = os.path.join(o.outdir, sample)
		if not os.path.exists(path_newbam): safe_mkdir(path_newbam)

		##call subsample function
		#unique
		subsampled_unique_bam = os.path.join(o.outdir, sample, sample+"_unique.bam")
		subsample_by_readlength(unique_bam, o.readlen_cutoff, subsampled_unique_bam)
		#multi
		subsampled_multi_bam = os.path.join(o.outdir, sample, sample+"_multi.bam")
		subsample_by_readlength(multi_bam, o.readlen_cutoff, subsampled_multi_bam)














