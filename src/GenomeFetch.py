#!/usr/bin/python

"""
GenomeFetch:

    A framework for fast fetching of regions from chromosomes
    as fasta files.
    
"""

_version = '1.3'

help_string = """
    GenomeFetch-%s

    USAGE:

      -o     : species (hsa, mmu), (default=hsa)
      -a     : assembly (e.g. hg19, mm9)
      -c     : chromosome (e.g. chr1, chrM, chr2_random)
      -f     : from coordinate (1-based, both ends inclusive)
      -t     : to coordinate (1-based, both ends inclusive)
      -s     : strand (e.g. 1 or -1)
      -d     : genome folder (e.g. /mnt/crick/sandberglab/genomes/hg19), instead of species and assembly
      \n""" % _version

import os, sys
from optparse import OptionParser

rcMapRNA = {'A':'U', 'C': 'G', 'G': 'C', 'U': 'A', 'T': 'A', 'N': 'N',
            'a':'u', 'c': 'g', 'g': 'c', 'u': 'a', 't': 'a', 'n': 'n'}
rcMapDNA = {'A':'T', 'C': 'G', 'G': 'C', 'U': 'A', 'T': 'A', 'N': 'N',
            'a':'t', 'c': 'g', 'g': 'c', 'u': 'a', 't': 'a', 'n': 'n'}

def revComp(seq, useRNA=False):
    return ''.join( [useRNA and rcMapRNA.get(b,'N') or rcMapDNA.get(b, 'N') for b in reversed(seq)])

class GenomeFetch:
    """Main Class For Fetching Genome Sequences.

    Administrator needs to set the correct path to directory
    containing .bases files in self.directories dictionary.

    cmdline interface usage:
    default organism (hsa) and assembly (hg17):
    GenomeFetch.py -c chr3 -f 100000 -t 1000010 -s 1
    specific organism and assembly:
    GenomeFetch.py -o mmu -a mm6 -c chr3 -f 100000 -t 1000010 -s 1

    module interface:
    >> import GenomeFetch
    >> gf = GenomeFetch.GenomeFetch() # dafualt specie and assmembly 
    >> gf.get_seq_from_to('chr12',100000,100010)
    AGCAGGTGCCT
    """
    specieslist = ['hsa','mmu','rn','gga','pto','cfa', 'ce','xenTro']
    
    defaultassembly = {"hsa":"hg19",
                       "mmu":"mm9",
                       "rn" :"rn4",
                       "gga":"galGal3",
                       "ce" :"ce6"
                       }

    defaultspecies = "hsa"
    defaultstrand = 1
    
    fasta_directories = {"hsa": {"hg19" : "/mnt/crick/sandberglab/genomes/hg19/"},
                         "mmu": {"mm8": '/mnt/crick/sandberglab/genomes/mm8_norandom/',
                                 "mm9"  : '/mnt/crick/sandberglab/genomes/mm9_norandom/',
                                 'cast_mm9': '/mnt/crick/danielr/casthybrid/CASTEi_genome/CAST_mm9/',
                                 'cast_mm9_ooref': '/mnt/crick/danielr/casthybrid/oocyte_ref/chr/',
                                 'cast_mm9_mgisnp': '/mnt/crick/danielr/casthybrid/mgi_snp/genomeCASTEi/'},
                         "rn" : {"rn4"  : '/mnt/crick/sandberglab/genomes/rn4/'},
                         "ce" : {"ce6"  : '/mnt/crick/sandberglab/genomes/ce6/'},
                         "cfa": {"canFam2" : '/mnt/crick/sandberglab/genomes/canFam2/'}
                         
                         }
    

    


    species = None
    assembly = None
    fileSeq = None
    openChromosome = None
    fasta_directory = ''
    
    def __init__(self, species=None,assembly=None, genomedir=None):
        '''
        Provide either species (e.g. hsa) and assembly (e.g. hg19), or genomedir (e.g. /mnt/crick/sandberglab/genomes/hg19)
        '''
        if genomedir is not None:
            # ignore species and assembly, just use the provided files
            self.fasta_directory = genomedir
        else:
        
            # set species
            if species is not None:
                if not species in self.specieslist:
                    sys.stderr.write("WARNING: selected species (%s) is not available\n" % species)
                    sys.stderr.write("AVAILABLE SPECIES:\n %s" % "\n".join([k for k in self.specieslist]))
                    return None
                self.species = species
            elif assembly is None:
                self.species = self.defaultspecies
            elif assembly in self.fasta_directories[self.defaultspecies]:
                self.species = self.defaultspecies
            else:
                # detect species from assembly
                species_from_assembly = dict((a,s) for s,a_arr in self.fasta_directories.items() for a in a_arr)
                self.species = species_from_assembly[assembly]
        
            # set assembly
            if assembly is None:
                self.assembly = self.defaultassembly[self.species]
            else:
                self.assembly = assembly
                if not self.fasta_directories[self.species].has_key(self.assembly):
                    sys.stderr.write("WARNING: selected assembly is not available")
                    return None
            
            self.fasta_directory = self.fasta_directories[self.species][self.assembly]

        self.__get_file_sizes()

    def __get_file_sizes(self, useFa=True):
        """Internal function that determines the chromosome length of all
        .bases files in a directory."""
        # reads in and stores the file lengths of the files in directory
        # -> the chromosome lengths
        self.files = os.listdir(self.fasta_directory)
        self.lengths = {}
        for file in self.files:
            if useFa:
                self.lengths[file] = os.path.getsize("%s%s" % \
                                                     (self.fasta_directory,\
                                                      file))

    def get_seq_from_to(self, chromosome,coordFrom, coordTo,strand=None, fLeaveOpen=False, Fa=True):
        """Main Function for fetching a genome sequence.
        Input: chromosome as chr4,
               coordFrom as integer, coordTo as integer
               strand as 1 or -1
               fLeaveOpen flag bool, indicates whether the file handle should
               be left open for repeatedly fetching from the same chromosome.
               """
        
        
        # chech chromosome input       
        self.chromosome = chromosome
        if not (os.path.exists('%s%s.fa' % (self.fasta_directory, chromosome))):
            print('WARNING: sequence  "%s%s.fa" does not exist.'% \
                  (self.fasta_directory, chromosome))
            return None

        # check strand input
        if strand is None:
            strand = self.defaultstrand
        if not isinstance(strand, int):
            mStrStrand = {'PLUS': 1, '+': 1, '1': 1, 'MINUS': -1, '-': -1, '-1': -1}
            if isinstance(strand, str):
                strand = mStrStrand[strand.upper()]
            else:
                strand = self.defaultstrand

        # chech coordinate input
        if not isinstance(coordFrom, int):
            coordFrom = int(coordFrom)
        if not isinstance(coordTo, int):
            coordTo = int(coordTo)
                
        # open file
        if not self.fileSeq or not self.openChromosome == chromosome:
            if self.fileSeq:
                self.fileSeq.close()
            self.fileSeq = open('%s%s.fa'% (self.fasta_directory,
                                            chromosome), 'r') # r+


        # get chromosome length
        ChromosomeLength =  self.lengths["%s.fa" % chromosome]
   

        # check that the requested region make sense
        if coordFrom > coordTo and strand == 1:
            sys.stderr.write('WARNING: coordTo is smaller than coordFrom and strand = 1(from %d, to %d)\n'%(\
                coordFrom, coordTo))
            return None
            
        # reverse coordinates if coordFrom > coordTo and strand is either -1 or not specified
        if coordFrom > coordTo and strand == - 1 or strand == None:
            tmpCoordFrom = coordFrom
            coordFrom = coordTo
            coordTo = tmpCoordFrom
            strand = -1
            
        # check that the requested region is within chromosome boundaries
        if coordTo > ChromosomeLength:
            sys.stderr.write('WARNING: tried to read past end of sequence (to %d, len %d)\n'%(\
                coordTo, ChromosomeLength))
            coordTo = ChromosomeLength
            
        if coordFrom > ChromosomeLength:
            sys.stderr.write('WARNING: tried to start reading past end of sequence (from  %d, len %d)\n'%(\
                coordTo, ChromosomeLength))
            return None
        
        if coordFrom < 0:
            sys.stderr.write('WARNING: tried to read past beginning of sequence (read from %d, start at 1)\n'%(coordFrom))
            coordFrom = 1

        # READ SEQUENCE FROM FILE
        # LINE LENGTH in nucleotides = 50 
        #
        baseOffset = 2 + len(chromosome) + (coordFrom-1) / 50 # to the start of sequence line and for all return strokes
        retbeforeEnd = (coordTo-1)/50 - (coordFrom-1) / 50
        ofsBase = baseOffset + coordFrom - 1
        self.fileSeq.seek(ofsBase, 0)
        rawSeq = self.fileSeq.read( abs( coordFrom - coordTo) + 1 + retbeforeEnd ).replace("\n","").upper()

        # reverse complement if desired
        if strand == -1:
            rawSeq = revComp(rawSeq)

        # a sanity check for a successful fetch
        assert('!' not in rawSeq)

        # if file will be read again, leave it open
        if not fLeaveOpen:
            self.fileSeq.close()
        else:
            self.openChromosome = chromosome

        # return sequence
        return rawSeq

def help():
    sys.stdout.write(help_string)


def DownloadAsembly(specie, assembly, toPath):
    """Downloads chromosome fasta files for specie and assmembly to the specified path."""
    import subprocess
    cmd = "wget http://hgdownload.cse.ucsc.edu/goldenPath/%s/bigZips/chromFa.tar.gz /tmp/infile.tar.gz" % assembly
    subprocess.call(cmd, shell=True)

    # unpack into correct path
    



if __name__=='__main__':
    opts = OptionParser()
    opts.add_option('-o','--species', dest='species')
    opts.add_option('-a','--assembly', dest='assembly') 
    opts.add_option('-c','--chromosome', dest='chromosome')    
    opts.add_option('-s','--seq-strand', dest='strand')    
    opts.add_option('-f','--seq-from', dest='coordFrom')        
    opts.add_option('-t','--seq-to',  dest='coordTo')    
    opts.add_option('-d','--genomedir', dest='genomedir')

    (options, args) = opts.parse_args()
    if options.chromosome is None or options.coordFrom is None or options.coordTo is None:
        help()
    else:
        gf = GenomeFetch(species=options.species, assembly=options.assembly, genomedir=options.genomedir)
        print gf.get_seq_from_to(options.chromosome,
                                 int(options.coordFrom),
                                 int(options.coordTo),
                                 strand=options.strand)
    
