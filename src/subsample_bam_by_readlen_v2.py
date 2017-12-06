import os
import argparse
import subprocess

def safe_mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)
        os.chmod(path, 0o774)

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
		bam = os.path.join(o.stardir, sample, "%s.bam" %sample)
		path_newbam = os.path.join(o.outdir, sample)
		if not os.path.exists(path_newbam): safe_mkdir(path_newbam)

		##call subsample function
		#unique
		subsampled_bam = os.path.join(o.outdir, sample, sample+".bam")
		subsample_by_readlength(bam, o.readlen_cutoff, subsampled_bam)














