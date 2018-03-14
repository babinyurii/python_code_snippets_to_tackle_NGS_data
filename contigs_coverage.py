#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 14:39:14 2018

@author: yuriy
"""

from Bio import SeqIO
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt
import seaborn as sns
import re
from math import log2
sns.set()



def contigs_cover_spades(file_name_or_path):
    """takes fasta file (spades tool output) 
    name or its path as an argument. 
    creates two plots:
    1. distribution of GC content in contigs
    2. GC content vs log2 coverage depth """
    
    container = []
    
    # using seq_record obj to access its attributes: 
    # id - name, seq - nucleotides
    for seq_record in SeqIO.parse(file_name_or_path, "fasta"):
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
    fig1 = plt.figure()
    plt.scatter(gc, cov_log2, s=5)
    plt.xlabel("GC content, %")
    plt.ylabel("log2 coverage depth")
    plt.title("contigs' coverage by reads")
    plt.savefig("GC_content_vs_contigs_coverage", format="jpeg")


# function call
contigs_cover_spades("/home/yuriy/sal/contigs_003.fasta")




    
    

    
    



