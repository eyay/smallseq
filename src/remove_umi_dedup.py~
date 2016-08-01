from __future__ import division
import os
import dr_tools, sys, argparse
import subprocess
from joblib import Parallel, delayed

def safe_mkdir(path):
	if not os.path.exists(path):
		os.mkdir(path)
		os.chmod(path, 0o774)

def remove_umi_dub(in_bam, out_bam, out_bed):
	cmd = 'umitools rmdup %s %s > %s' %(in_bam, out_bam, out_bed)
	subprocess.call(cmd, shell = True)

def dedup(in_bam, out_bam, logfile, prefix):
	cmd = 'python src/public/UMI-tools/dedup_umi.py --method %s --output-stats %s --further-stats --read-length -I %s -S %s -L %s' %(o.method, prefix, in_bam, out_bam, logfile)
	subprocess.call(cmd, shell = True)

def parallel_process(sample):
	###prepare input files
	#bam files
	unique_bam = os.path.join(o.stardir, sample, "%s_unique.bam" %sample)
	multi_bam = os.path.join(o.stardir, sample, "%s_multi.bam" %sample)

	###prepare output files
	path_newbam = os.path.join(o.outdir, sample)
	if not os.path.exists(path_newbam): safe_mkdir(path_newbam)
	#bam files
	UMIdubRemoved_unique_bam = os.path.join(o.outdir, sample, sample+"_unique.bam")
	UMIdubRemoved_multi_bam = os.path.join(o.outdir, sample, sample+"_multi.bam")
	#log files
	path_logfiles = os.path.join(o.outdir, sample, 'logs')
	if not os.path.exists(path_logfiles): safe_mkdir(path_logfiles)
	unique_log = os.path.join(path_logfiles, "unique_%s.log" %(o.method))
	multi_log = os.path.join(path_logfiles, "multi_%s.log" %(o.method))
	#prefix path
	path_stats = os.path.join(o.outdir, sample, 'stats')
	if not os.path.exists(path_stats): safe_mkdir(path_stats)
	unique_prefix = os.path.join(path_stats, "unique_%s" %(o.method))
	multi_prefix = os.path.join(path_stats, "multi_%s" %(o.method))
	#call function
	dedup(unique_bam, UMIdubRemoved_unique_bam, unique_log, unique_prefix)
	dedup(multi_bam, UMIdubRemoved_multi_bam, multi_log, multi_prefix)


if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--stardir', required=True)
	parser.add_argument('-o', '--outdir', required=True)
	parser.add_argument('-m', '--method', default="adjacency")
	parser.add_argument('-p', '--numCPU', default=20, type=int)
	o = parser.parse_args()

	if not os.path.exists(o.outdir): safe_mkdir(o.outdir)

	sample_names = os.listdir(o.stardir)
	Parallel(n_jobs=o.numCPU)(delayed(parallel_process)(sample) for sample in sample_names)









