#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 14:13:54 2018
@author: yuriy
"""

from Bio import SeqIO


def cat_contigs(file_or_path, seq_id="long_seq", output ="long_seq.fasta"):
    """creating long sequence 
    from multifasta file containing contigs. 
    usefule to blast the resulting seq
    args: 1. file name or path to file
          2. name of the seq (after > in fasta file)
          3. output - output file name """
    
    long_seq = ""
    for record in SeqIO.parse(file_or_path, "fasta"):
        long_seq += record
   
    long_seq.id = seq_id
    SeqIO.write(long_seq, output, "fasta")


cat_contigs("/home/yuriy/sal/contigs_003.fasta", seq_id="my_func_seq_output", output="my_func_out.fasta")


    
    
    
    
    



    


    


    