# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 11:33:22 2019

@author: babin
"""
import os
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio import Entrez
import datetime
import seaborn as sns
from time import sleep, time
sns.set()

def _get_current_time():
    """just returns time stamp
    """
    time_stamp = datetime.datetime.fromtimestamp(
        time()).strftime('%Y-%m-%d %H:%M:%S')
    return time_stamp

def _genbank_loader(f_obj, seq_id, rettype):
    handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype=rettype, retmode="text")
    fetched = handle.read()
    f_obj.write(fetched)


def fetch_seq(ids, seq_format="fasta", sep=False):
    """downloads sequences from nucleotide database
    by id nums and saves them in genbank format
    ----------
    ids : str or list of str
        sequence genbank id or list of ids
    seq_format : str
        gb - genbank files
        fasta (by default) - fasta files
    sep : bool
        False - download bunch of sequences as one file
        True - donwload bunch of sequences as separate files
    """
    Entrez.email = "babin.yurii@gmail.com"
    count = 0
    if type(ids) == str:
        with open("downloaded_" + ids + "." + seq_format, "w") as f_obj:
            _genbank_loader(f_obj, ids, seq_format)
            print("a sequence " + ids + " was downloaded")
    elif type(ids) == list:
        if sep:
            for i in ids: 
                with open("downloaded_" + i + "." + seq_format, "w") as f_obj:
                    _genbank_loader(f_obj, i, seq_format)
                    count += 1
                    sleep(0.5)
            print("a total of %s sequences were downloaded" %count)
        else:
            time_stamp = _get_current_time()
            days, day_time = time_stamp.split(" ")
            day_time = day_time.split(":")
            day_time = "_".join(day_time)
            time_stamp = days + "_time-" + day_time
            with open("downloaded_bunch_" + time_stamp + "." + seq_format, "w") as f_obj:
                for i in ids:
                    _genbank_loader(f_obj, i, seq_format)
                    count += 1
                    sleep(0.5)
                print("a total of %s sequences were downloaded" %count)
    else:
        print("invalid ids parameter type")


def _get_id_length_gc(file):
    """calculates length and GC content of the sequences
    Parameters
    ----------
    file : input fasta file
    
    Returns
    ----------
    num_records : number of sequences in the input fasta file
    ids_len_and_gc : list of tuples, each tuple contains sequence id, length and GC content
    """
    ids_len_and_gc = []
    records = SeqIO.parse(file, "fasta")
    num_records = 0
    for rec in records:
        ids_len_and_gc.append((rec.id, len(rec.seq), GC(rec.seq)))
        num_records += 1
    return num_records, ids_len_and_gc
        
        
  
def _show_fasta_info(file, num_records, ids_len_and_gc):
    """prints out the result of fasta_info function
    Parameters
    ----------
    file : input fasta file
    num_records : number of sequences in the input file
    ids_len_and_gc : list of tuples containing sequence id, length and GC content
    """
    print("file '{0}' contains {1} sequences".format(file, num_records))
    print("", "sequence id", "length", "GC%", sep="\t")
    for counter, value in enumerate(ids_len_and_gc, 1):
        print(counter, value[0], value[1], round(value[2], 2), sep="\t")
        print("------------------------------------")
        
def fasta_info(path_to_fasta=False):
    """prints out information about fasta files:
    number of sequences in the file, sequence id numbers,
    lengths of sequences and GC content
    
    without arguments takes as an input
    all fasta files in the current dir
    
    Parameters
    ----------
    path_to_fasta : str or list
        path to input file, or list of paths
    """
    fasta_extensions = ["fa", "fas", "fasta"]
    
    if type(path_to_fasta) == str:
        num_records, len_and_gc = _get_id_length_gc(path_to_fasta)
        _show_fasta_info(path_to_fasta, num_records, len_and_gc)
        
    elif type(path_to_fasta) == list:
        for path in path_to_fasta:
            num_records, len_and_gc = _get_id_length_gc(path)
            _show_fasta_info(path, num_records, len_and_gc)  
    else:
        current_dir_content = os.listdir()
        for f in current_dir_content:
            if f.rsplit(".", 1)[-1] in fasta_extensions:
                num_records, ids_len_and_gc = _get_id_length_gc(f)
                _show_fasta_info(f, num_records, ids_len_and_gc)





























      
        
        
        
        
        
        
        
