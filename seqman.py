#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 30 15:13:31 2018

@author: yuriy
"""

from Bio import SeqIO
from Bio.SeqUtils import GC
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import re
from math import log2
sns.set()



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





def bowtie2_run(bowtie_build, ref_seq, prefix, bowtie_align, mapping_name, 
                fastq_file, output, summary):
    """runs bowtie2
    
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
     
        


def fastqc_run(path_to_fastqc, path_to_file, file_log=False):
    """runs fastqc  
    
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


def prefetch_run(path_to_prefetch, sra_ids):
    """
    launch prefetch utility
    
    args:
    path to prefetch util
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
            
            


def fastq_dump_run(path_to_sra, out_dir=False):
    """
    run fastq-dump utility
    
    args:
    path to sra util
    path to output folder
    """
    
    if out_dir:
        
        with open(path_to_sra) as f_obj:
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
                
         
        
        
        
        
