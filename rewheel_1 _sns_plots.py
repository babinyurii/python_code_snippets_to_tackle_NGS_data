#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 11:09:32 2018

@author: yuriy
"""
from Bio import SeqIO
from statistics import mean
from time import time
import matplotlib.pyplot as plt
from Bio.SeqUtils import GC
import seaborn as sns
sns.set()


# for seaborn distributions:
# sns.kdeplot(data, shade=True)
# but the better is sns.distplot(data)


def count_reads(input_file, file_ext="fastq"):
    
    counter = 0
    for record in SeqIO.parse(input_file, file_ext):
        counter += 1
    return counter

# the same as count_reads but generator function
def count_reads_gen(input_file, file_ext="fastq"):
    
    for record in SeqIO.parse(input_file, file_ext):
        yield record


def longest_shortest_read(input_file, file_ext="fastq"):
    lens = [len(record) for record in SeqIO.parse(input_file, file_ext)]
    seq_gc = [GC(record.seq) for record in SeqIO.parse(input_file, file_ext)]  # from here your can take mean GC value
    phred_per_seq = [mean(record.letter_annotations["phred_quality"]) \
                       for record in SeqIO.parse(input_file, file_ext)]
    
    lens_dist(lens)
    #return (min(lens), max(lens), mean(lens))
    gc_dist(seq_gc)
    phred_dist(phred_per_seq)
    

def lens_dist(lens_data):
    fig_dist = plt.figure()
    sns.distplot(lens_data).set_title("length distribution")
    plt.show()


def gc_dist(seqs_data):
    fig_gc_dist = plt.figure()
    sns.distplot(seqs_data).set_title("GC distribution")
    plt.show()
    
def phred_dist(phred_data):
    fig_phred = plt.figure()
    sns.distplot(phred_data).set_title("phred quality distribution per read")
    plt.show()



# call
print(longest_shortest_read("/home/yuriy/sal/sal_data/sal_raw_003.fastq"))





























