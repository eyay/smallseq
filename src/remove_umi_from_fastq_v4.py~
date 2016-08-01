from __future__ import division
import os
import dr_tools, sys, argparse
import subprocess
from joblib import Parallel, delayed

"""
parallelized version
this v4 is parallelized and uses CGAT umi removal tool
"""

def safe_mkdir(path):
	if not os.path.exists(path):
		os.mkdir(path)
		os.chmod(path, 0o774)

def trim_umi(in_fq, out_fq, umi_seq, logfile):
	cmd = 'python src/public/UMI-tools/extract_umi.py --bc-pattern=%s -L %s -I %s -S %s' %(umi_seq, logfile, in_fq, out_fq)
	subprocess.call(cmd, shell = True)

def parallel_process(sample):
	fq_file_list = os.listdir(os.path.join(rawdatadir, sample))
	for afile in fq_file_list:
		raw_fq = os.path.join(rawdatadir, sample, afile)
	##prepare output files
	path_trimmed_fq = os.path.join(outdir, sample)
	if not os.path.exists(path_trimmed_fq): safe_mkdir(path_trimmed_fq)

	fq_name = raw_fq.split("/")[-1].split(".",1)[0]
	trimmed_fq = os.path.join(outdir, sample, "%s_umiTrim.fq" %fq_name)
	logfile = os.path.join(outdir, sample, "extract.log")

	##call subsample function
	trim_umi(raw_fq, trimmed_fq, umi_seq, logfile)


if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--rawdatadir', required=True)
	parser.add_argument('-o', '--outdir', required=True)
	parser.add_argument('-u', '--umi_seq', default="HHHHHHCA", type=str)
	parser.add_argument('-p', '--numCPU', default=20, type=int)
	o = parser.parse_args()

	if not os.path.exists(o.outdir): safe_mkdir(o.outdir)
	rawdatadir = o.rawdatadir
	outdir = o.outdir
	umi_seq = o.umi_seq

	sample_names = os.listdir(o.rawdatadir)
	Parallel(n_jobs=o.numCPU)(delayed(parallel_process)(sample) for sample in sample_names)











