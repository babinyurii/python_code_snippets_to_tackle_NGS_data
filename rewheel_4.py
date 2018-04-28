#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 28 16:48:32 2018

@author: yuriy
"""


from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()




contig_lens = []
counter = 0
for record in SeqIO.parse("/home/yuriy/sal/sal_data/contigs_003.fasta", "fasta"):
    counter += 1
    contig_lens.append(len(record))



total_len = sum(contig_lens)
max_contig = max(contig_lens)
min_contig = min(contig_lens)   


fig_hist = plt.figure()
sns.distplot(contig_lens, norm_hist=False)
plt.show()



# the same but with raw counts
fig_hist_1 = plt.figure()
sns.distplot(contig_lens, kde=False, hist=True)



