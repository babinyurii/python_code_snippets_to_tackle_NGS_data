                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                

from Bio import SeqIO
from statistics import mean
from time import time
import matplotlib.pyplot as plt
from Bio.SeqUtils import GC
import seaborn as sns
sns.set()

from statistics import mean, stdev


def phred_per_base(input_file, file_ext="fastq"):
    
    lens = [len(record) for record in SeqIO.parse(input_file, file_ext)]
    max_len = max(lens)
    d_phred = {k:() for k in range(1, max_len + 1)}
    
    for rec in SeqIO.parse(input_file, file_ext):
        ph = rec.letter_annotations["phred_quality"]
        counter = 1
        for i in ph:
            d_phred[counter] += (i,) # under the key increment the tuple
            counter += 1
    return d_phred


d = phred_per_base("/home/yuriy/sal/sal_data/sal_raw_003.fastq")

phred_mean = []
phred_stdev = []

for key, value in d.items():
    phred_mean.append(mean(value))
    phred_stdev.append(stdev(value + (0,)))
    
    
fig_bar_phred = plt.figure(figsize=(10, 5))
plt.bar(list(range(len(phred_mean))), phred_mean)
plt.show()

fig_scatter_phred = plt.figure(figsize=(10,8))
plt.scatter(list(range(len(phred_mean))), phred_mean)
plt.show()

fig_plot_phred = plt.figure(figsize=(10,8))
plt.plot(phred_mean)
plt.axhline(20, c="red")
plt.show()