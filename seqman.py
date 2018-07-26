#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 30 15:13:31 2018
@author: yuriy
"""

from Bio import SeqIO
from Bio import AlignIO
from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.SeqUtils import GC
import glob
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import re
from math import log2
from operator import itemgetter
from random import randint
from time import sleep, time
sns.set()


def blast_fasta(path_to, e_thresh=0.1, hits_to_return=10):  
    """
    takes fasta file and blasts all the record found in it
    writes results into the txt file
    args:
    path_to - path to fasta file
    e_thresh - e-value cut off
    hits_to_return - a number of hits returned
    """
    fasta_to_blast = SeqIO.parse(path_to, "fasta")
    
    with open("blast_results_seqman.txt", "w") as f_obj:
        start = time()
        for rec in fasta_to_blast:
            f_obj.write("\nQUERY: " + rec.id + "\n")
        
            result_handle = NCBIWWW.qblast("blastn", "nt",  rec.seq, hitlist_size=hits_to_return)
            blast_record = NCBIXML.read(result_handle)
        
            for alignt in  blast_record.alignments:
                for hsp in alignt.hsps:
                    if hsp.expect < e_thresh:
                        f_obj.write("---------hit--------\n")
                        f_obj.write("sequence: " + alignt.title + "\n")
                        f_obj.write("length: " + str(alignt.length) + "\n")
                        f_obj.write("e value: " + str(hsp.expect) + "\n")
               
            end = time()
            print(rec.id + " blast query was finished in {0} minutes {1} seconds".format((end - start) // 60, int((end - start) % 60)))
            
            
def fasta_info(path_to):
    """
    returns information about fasta file
    ids and length of sequences
    ----------
    note: 
    if fasta in the current location, use .\file.fasta
    or ./file.fasta as path_to arg
    ----------
    """
    container = []
    records = SeqIO.parse(path_to, "fasta")
    
    for rec in records:
        container.append((rec.id, len(rec.seq)))
    print("no", "id", "length")
    print("------------------")
    for counter, value in enumerate(container, 1):
        print(counter, value[0], value[1])
    #return container


def fetch_by_id(ids):
    """
    downloads sequences from nucleotide database
    by id nums
    creates 2 files: .gb and .fasta
    ids - list of ids to download
    returns no value
    """
    Entrez.email = "babin.yurii@gmail.com"
    
    with open("downloaded_seq.gb" , "w") as f_obj:
        for i in ids: 
            sleep(0.5)
            handle = Entrez.efetch(db="nucleotide", id=i, rettype="gb", retmode="text")
            fetched = handle.read()
            f_obj.write(fetched)
    
    count = 0
    with open("downloaded_seq.gb", "r") as input_handle:
        with open("downloaded_seq.fasta", "w") as output_handle:
            seqs = SeqIO.parse(input_handle, "genbank")
            for seq in seqs:
                SeqIO.write(seq, output_handle, "fasta")
                count += 1
    print("a total of %s sequences were downloaded" %count)
         
    
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
    splits multifasta into several fastas
    saves the results in the original fasta dir 
    ----------
    note: 
    if multifasta in the current location, use .\multifasta.fasta
    or ./multifasta.fasta as path_to_file arg
    ----------
    """
    path = path_to_file.rsplit("\\", 1)[0]
    for record in SeqIO.parse(path_to_file, "fasta"):        
        SeqIO.write(record, path + "\\"  + record.id + ".fasta", "fasta")


def merge_fasta(folder_path):
    """
    merges fastas with single or multi records
    into a multifasta file
    takes path to folder as an arg
    -------
    note: if fastas are in the current location,
    use ".\\" or "./" as folder_path arg
    -------
    """  
    fastas = glob.glob(folder_path + "*.fasta")
    with open(folder_path + "merged_multifasta.fasta", "w") as f_obj:
        for el in fastas:
            try:
                rec = SeqIO.read(el, "fasta")
                SeqIO.write(rec, f_obj, "fasta")
            except ValueError as e:
                print("in the %s file more than one record exist" %el)
                multi_recs = SeqIO.parse(el, "fasta")
                for r in multi_recs:
                    SeqIO.write(r, f_obj, "fasta")
    
   
    
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
    





        
        
