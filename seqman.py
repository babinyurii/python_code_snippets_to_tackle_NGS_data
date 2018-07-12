#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 30 15:13:31 2018
@author: yuriy
"""

from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqUtils import GC
import glob
import subprocess
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import re
from math import log2
from operator import itemgetter
from random import randint
sns.set()

def genseq(seq_len):
    """
    the most simplest 
    random sequence generator
    seq_len parameter is a desired
    sequence length
    """
    nucs = ["A", "T", "G", "C"]
    seq = ""
    for i in range(0, seq_len):
        seq += nucs[randint(0,3)]
    return seq


def cat_contigs(path_to, seq_id="long_seq", out_f ="long_seq.fasta"):
    """creating long sequence 
    from multifasta file containing contigs. 
    usefule to blast the resulting seq.
    args: 1. path to file
          2. name of the seq (">name_seq" in fasta)
          3. output file path and name """
    
    cat_seq = ""
    for record in SeqIO.parse(path_to, "fasta"):
        cat_seq += record
   
    cat_seq.id = seq_id
    SeqIO.write(cat_seq, out_f, "fasta")


def split_fasta(path_to_file):
    """
    splits multifasta and writes resulting sequences 
    into several fasta files
    saves the results in the dir of the original multifasta file 
    """
    path = path_to_file.rsplit("\\", 1)[0]
    for record in SeqIO.parse(path_to_file, "fasta"):        
        SeqIO.write(record, path + "\\"  + record.id + ".fasta", "fasta")


def merge_fasta(folder_path):
    """
    merges several fasta with single record in each
    into a multifasta
    takes path to folder as an argument
    all the .fasta in this folder will be merged
    if there's multirecord file, throws: 
    "ValueError: More than one record found in handle"
    """
    
    fastas = glob.glob(folder_path + "*.fasta")
    with open(folder_path + "merged_multifasta.fasta", "a") as f_obj:
        for el in fastas:
            rec = SeqIO.read(el, "fasta")
            SeqIO.write(rec, f_obj, "fasta")
    
   
    
def locate_gaps(path_to, gaps=4, exact=True):
    """return the position and column in an alignment
    with the give number of gaps and less 
    
    path_to file
    gaps - optional number of gaps to find
    
    now reads only fasta
    """
    
    alignment = AlignIO.read(path_to, "fasta")
    num_cols = len(alignment[0])
    container = []
    

    if exact:
        for i in range(num_cols):
            col = alignment[ : , i]
            num_gaps = col.count("-")
            if num_gaps == gaps:
                container.append((i + 1, col))
    
    else:  
        for i in range(num_cols):
            col = alignment[ : , i]
            num_gaps = col.count("-")
            if 1 <= num_gaps <= gaps:
                container.append((i + 1, col))
     
    
    return container


def locate_ambig(path_to, line=False):
    """
    return the position and column containing 
    ambiguous bases in list of tuples
    
    args:
    path to fasta alignment, 
    'if line=int', find ambiguous nucleotides 
    at a certain seq in an alignment
    'if line=False, finds all the ambiguous nucs
    """
    ambig_dna = ["Y", "R", "W", "S", "K", "M", "D", "V", "H", "B", "X", "N"]
    alignment = AlignIO.read(path_to, "fasta")
    container = []

    if line:
        num_cols = len(alignment[0])    
    
        for i in range(num_cols):
            col = alignment[line, i]
            for nuc in ambig_dna:
                if nuc in col:
                    container.append((i + 1, col))
                    break
    else:
        num_cols = len(alignment[0])
    
        for i in range(num_cols):
            col = alignment[ : , i]
            for nuc in ambig_dna:
                if nuc in col:
                    container.append((i + 1, col))
                    break
    
    return container
    
 
def filter_gc_depth(input_file, gc_low, gc_high, cov_low, cov_high):
    """filtering contigs by GC content and coverage depth
    creates fasta file containing resulted contigs. 
    for input spades file
    
    input file is fasta
    gc and cov are integers
    """
    fraction = []
    gc = []
    #coverage = []
    reg = re.compile(r'cov_(.*)')
    
    for seq_record in SeqIO.parse(input_file, "fasta"):
        
        gc = GC(seq_record.seq)
        cov_raw = reg.search(seq_record.id)
        cov = cov_raw.group()
        cov = log2(float(cov[4:]))
        
        if gc_low <= gc <= gc_high and cov_low <= cov <= cov_high:
            fraction.append(seq_record)
        
    SeqIO.write(fraction, "fraction.fasta", "fasta")
    
    
def contigs_cover_spades(path_to_file):
    """takes fasta file (spades tool output) 
    name or its path as an argument. 
    creates two plots:
    1. distribution of GC content in contigs
    2. GC content vs log2 coverage depth """
    
    container = []
    
    # using seq_record obj to access its attributes: 
    # id - name, seq - nucleotides
    for seq_record in SeqIO.parse(path_to_file, "fasta"):
        entry = ()
        entry = (seq_record.id, GC(seq_record.seq))
        container.append(entry)
        entry = ()
      
    gc = [x[1] for x in container]
    
    # fig1 = plt.figure()
    sns_dist = sns.distplot(gc, hist=False, 
                            kde_kws={"shade":True})
    sns_dist.set_title("GC_distribution")
    sns_dist.set_xlabel("GC content, %")
    fig = sns_dist.get_figure()
    fig.savefig("contigs_GC_distribution", format="jpeg")
    
    
    reg = re.compile(r'cov_(.*)')
    coverage = []
    for el in container:
        cov_raw = reg.search(el[0]) 
        cov = cov_raw.group()
        cov = cov[4:]
        coverage.append(float(cov))
    
    cov_log2 = [log2(x) for x in coverage]
    fig1 = plt.figure(figsize=(10, 8))
    plt.scatter(gc, cov_log2, s=5)
    plt.xlabel("GC content, %")
    plt.ylabel("log2 coverage depth")
    plt.title("contigs' coverage by reads")
    plt.savefig("GC_content_vs_contigs_coverage", format="jpeg")
        


def coverage_count(ugene_cov):
    """
    takes ugene export with the following parameters:
    format: 
        - per base (it derives .txt file, separator is tab)
        - export bases quantity 'yes'
    ----
    sorts bases count at each position
    useful to verify ambigous nucleotides
    writes the result into a .csv in the current location, 
    named 'coverage_depth_count.csv"
    returns no container
    """
    
    cov = pd.read_csv(ugene_cov, sep="\t")
    cov.index = cov.position
    #cov.drop(["position"], axis=1, inplace=True)
    
    nucleotides = cov.columns.tolist()[-4:]
    container = []

    for ind in cov.index:
        nucl_count = cov.loc[ind, "A" : "T"].tolist()
        count_at_pos = list(zip(nucleotides, nucl_count))
        count_at_pos.sort(key=itemgetter(1), reverse=True)
        # creating more nicely looking thing for resulted series obj:
        
        s = ""
        for el in count_at_pos:
            s += el[0] + ": " + str(el[1]) + ",     "
            
        container.append(s)
    
    cov["coverage_count"] = container  
    cov.to_csv("coverage_depth_count.csv")
    

def bowtie2_run(bowtie_build, ref_seq, prefix, bowtie_align, mapping_name, 
                fastq_file, output, summary):
    """
    runs bowtie2
    
    args: 
    path to bowtie_build.py
    path to reference seq (fasta)
    path to and prefix file name (only id letters)
    name of the run (for log)
    input fastq file
    output .sam file name
    path and name.txt for summary
    """
    
    build_args = [bowtie_build, ref_seq, prefix]
    
    build_inst = subprocess.Popen(build_args, stderr=subprocess.PIPE)
    for line in build_inst.stderr:
        print(line)
        
        
    align_args = [bowtie_align, "-x", prefix, "-U", fastq_file, "-S", output]
    
    align_inst = subprocess.Popen(align_args, stderr=subprocess.PIPE)
    
    with open(summary, "a") as f_obj:
        f_obj.write(mapping_name + "\n")
        for line in align_inst.stderr:
            f_obj.write(str(line) + "\n")
        f_obj.write("***********end of run************")
        
        

def fastqc_run(path_to_fastqc, path_to_file, file_log=False):
    """
    runs fastqc  
    
    'path_to_fastqc' arg should be a file 
    name, path, or a list of names/paths 
    'path_to_file' is path to fastq file
    if file_log=True, stderr of fastqc
    will be written into a file
    """
    
    if type(path_to_file) == list:
        for name_path in path_to_file:
            fastqc_inst = subprocess.Popen([path_to_fastqc, name_path], 
                                           stderr=subprocess.PIPE)
            
            if file_log:
                with open("fastqc_log.txt", "a") as f_obj:
                    for line in fastqc_inst.stderr:
                         f_obj.write(str(line) + "\n")
                    f_obj.write("***************\n")
            
            else:
                for line in fastqc_inst.stderr:
                    print(line)
                
    else:
        fastqc_inst = subprocess.Popen([path_to_fastqc, path_to_file], 
                                       stderr=subprocess.PIPE)
        
        
        # fastqc_inst.kill() # to kill the process
        
        if file_log:
            with open("fastqc_log.txt", "a") as f_obj:
                for line in fastqc_inst.stderr:
                    f_obj.write(str(line) + "\n")
                f_obj.write("***************\n")
        else:
             for line in fastqc_inst.stderr:
                    print(line)



def prefetch_run(sra_ids):
    """
    launch prefetch utility
    
    args:
    path to sra ids in the txt file
    """
    
    with open(sra_ids) as f_obj:
        for line in f_obj:
            line = line.strip()
            prefetch_inst = subprocess.Popen(["/home/yuriy/tools/sratoolkit.2.9.0-centos_linux64/bin/prefetch", line],
                                            stderr=subprocess.PIPE)
            
            for line in prefetch_inst.stderr:
                print(line)
                print("**********************")
            
            


def fastq_dump_run(path_to, out_dir=False):
    """
    run fastq-dump utility
    path to the utility itself is in fastq_dump_inst, 
    change path if needed    

    args:
    path to file with accessions
    path to output folder
    """
    
    if out_dir:
        
        with open(path_to) as f_obj:
            for line in f_obj:
                line = line.strip()
                fastq_dump_inst = subprocess.Popen(["/home/yuriy/tools/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump",
                                                    "--outdir",
                                                    out_dir,
                                                    "--split-files",
                                                    line],
                                                    stderr=subprocess.PIPE)
                for line in fastq_dump_inst.stderr:
                    print(line)
    
    
    else:
        with open(path_to_sra) as f_obj:
            for line in f_obj:
                fastq_dump_inst = subprocess.Popen(["/home/yuriy/tools/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump",
                                                    "--split-files",
                                                    line],
                                                    stderr=subprocess.PIPE)
                for line in fastq_dump_inst.stderr:
                    print(line)
                
         
        
        
        
        
