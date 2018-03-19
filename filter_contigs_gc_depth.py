#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 12:08:00 2018
@author: yuriy
"""

from Bio import SeqIO
from Bio.SeqUtils import GC
import re
from math import log2


def filter_gc_depth(input_file, gc_low, gc_high, cov_low, cov_high):
    """filtering contigs by GC content and coverage depth
    creates fasta file containing resulted contigs
    
    input file is fasta
    gc and cov are integers
    """
    fraction = []
    gc = []
    coverage = []
    reg = re.compile(r'cov_(.*)')
    
    for seq_record in SeqIO.parse(input_file, "fasta"):
        
        gc = GC(seq_record.seq)
        cov_raw = reg.search(seq_record.id)
        cov = cov_raw.group()
        cov = log2(float(cov[4:]))
        
        if gc_low <= gc <= gc_high and cov_low <= cov <= cov_high:
            fraction.append(seq_record)
        
    SeqIO.write(fraction, "fraction.fasta", "fasta")
    
   
    
filter_gc_depth("contigs_004.fasta", 40, 60, 0, 4)


        