#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 28 15:33:54 2018

@author: yuriy
"""
from Bio import SeqIO
from statistics import mean
from time import time
import matplotlib.pyplot as plt
from Bio.SeqUtils import GC
import seaborn as sns
sns.set()

from statistics import mean, stdev



def nucleotide_per_base(input_file, file_ext="fastq"):
    
    lens = [len(record) for record in SeqIO.parse(input_file, file_ext)]
    max_len = max(lens)
    d_seq = {k:() for k in range(1, max_len + 1)}    # make frozen dict which has order!!!!
    
    for rec in SeqIO.parse(input_file, file_ext):
        s = rec.seq
        print(s)
        counter = 1
        for char in s:
            d_seq[counter] += (char,) # under the key increment the tuple
            counter += 1
    return d_seq

d = nucleotide_per_base("/home/yuriy/sal/sal_data/sal_raw_003_short_for_test.fastq")



gc_per_base = []
for k,v in d.items():
    gc_per_base.append(GC(v))
    
plt.plot(gc_per_base)



# generator GC per base
gc_per_base_gen = (GC(s) for k,s in d.items())


# generator, contains sequences
seqs_gen = (s for k,s in d.items())


