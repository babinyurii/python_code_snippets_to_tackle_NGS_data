# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 11:33:22 2019

@author: babin
"""
import os
import datetime
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from time import sleep, time
import matplotlib.pyplot as plt
import seaborn as sns
from math import log2
from http.client import IncompleteRead
from socket import gaierror
from urllib.error import HTTPError
import pandas as pd
sns.set()


def _get_current_time():
    time_stamp = datetime.datetime.fromtimestamp(
        time()).strftime('%Y-%m-%d %H:%M:%S')
    return time_stamp

def _load_from_genbank(f_obj, seq_id, rettype):
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
            _load_from_genbank(f_obj, ids, seq_format)
            print("a sequence " + ids + " was downloaded")
    elif type(ids) == list:
        if sep:
            for i in ids: 
                with open("downloaded_" + i + "." + seq_format, "w") as f_obj:
                    _load_from_genbank(f_obj, i, seq_format)
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
                    _load_from_genbank(f_obj, i, seq_format)
                    count += 1
                    sleep(0.5)
                print("a total of %s sequences were downloaded" %count)
    else:
        print("invalid ids parameter type")



def _fetch_blast_results(record, e_thresh, hits):
    result_handle = NCBIWWW.qblast("blastn", "nt",  record.seq, hitlist_size=hits)
    blast_record = NCBIXML.read(result_handle)
    blast_results_record = []
    for alignt in  blast_record.alignments:
        for hsp in alignt.hsps:
            if hsp.expect < e_thresh:
                blast_results_record.append([record.id, alignt.title, str(alignt.length), str(hsp.expect)])
    return blast_results_record
                


def blast_fasta(query, e_thresh=0.1, hits=1):
    """blast records from a fasta file
    writes results into the tab-delimited txt file
    
    Parameters:
    -----------
    query: str
        path to the input file
    e_thresh: float
        e-value blast threshold 
    hits: int
        a number of hits to return, 1 by default
    """
    fasta = SeqIO.parse(query, "fasta")
    blast_results_total = []
    
    for record in fasta:
        try:
            blast_results_record = _fetch_blast_results(record, e_thresh, hits)
            for res in blast_results_record:
                blast_results_total.append(res)
                
            time_stamp = _get_current_time()
            print(record.id, " blasted at: ", time_stamp)
        
        except IncompleteRead as e: 
            print("Network problem: ", e, "Second and final attempt is under way...")
            blast_results_record = _fetch_blast_results(record, e_thresh, hits)
            for res in blast_results_record:
                blast_results_total.append(res)
                
            time_stamp = _get_current_time()
            print(record.id, " blasted at: ", time_stamp)
            
        except gaierror as e:
            print("some other problem, 'gaierror': ", e)
            
        except HTTPError as e:
            print("urllib.error.HTTPError: ", e)
    
    df = pd.DataFrame(blast_results_total, columns=["record_id", "hit_name", "hit_length", "e_value"])
    df.to_csv("blast_results.csv", sep="\t")
    print("job done. the results are in {0}".format(os.path.abspath("blast_results.csv")))
            


def _get_id_length_gc(file):
    ids_len_and_gc = []
    records = SeqIO.parse(file, "fasta")
    num_records = 0
    for rec in records:
        ids_len_and_gc.append((rec.id, len(rec.seq), GC(rec.seq)))
        num_records += 1
    return num_records, ids_len_and_gc
        
        
  
def _show_fasta_info(file, num_records, ids_len_and_gc):
    print("file '{0}' contains {1} sequences".format(file, num_records))
    print("", "sequence id", "length", "GC%", sep="\t")
    for counter, value in enumerate(ids_len_and_gc, 1):
        print(counter, value[0], value[1], round(value[2], 2), sep="\t")
        print("------------------------------------")
        
        
def fasta_info(path_to=False):
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
    
    if type(path_to) == str:
        num_records, len_and_gc = _get_id_length_gc(path_to)
        _show_fasta_info(path_to, num_records, len_and_gc)
        
    elif type(path_to) == list:
        for path in path_to:
            num_records, len_and_gc = _get_id_length_gc(path)
            _show_fasta_info(path, num_records, len_and_gc)  
    else:
        current_dir_content = os.listdir()
        for f in current_dir_content:
            if f.rsplit(".", 1)[-1] in fasta_extensions:
                num_records, ids_len_and_gc = _get_id_length_gc(f)
                _show_fasta_info(f, num_records, ids_len_and_gc)


def split_fasta(path_to, path_out=False):
    """splits fasta file containing several
    sequences into the corresponding number of
    fasta files. 
    Parameters:
    ----------
    path to : str 
        path to the input file
    path_out : str
        path to output dir
    """
    if path_out:
        if not os.path.exists(path_out):
            os.mkdir(path_out)
        for record in SeqIO.parse(path_to, "fasta"):        
            SeqIO.write(record, path_out + record.id + ".fasta", "fasta")
        print("file {0} was splitted. the results are in the {1}".format(path_to, path_out))
    else:
        for record in SeqIO.parse(path_to, "fasta"):
            SeqIO.write(record, record.id + ".fasta", "fasta")
        print("file {0} was splitted. the results are in the {1}".format(path_to, os.getcwd()))


def _cat_fasta_records(file):
    cat_seq = ""
    for record in SeqIO.parse(file, "fasta"):
        cat_seq += record
    return cat_seq


def cat_fasta(path_to, fas_name="cat_seq.fasta", fas_id="cat_seq", fas_descr=""):
    """concatenates fasta sequences into one 
    long sequence, takes one multifasta 
    or several fasta files as an input
    Parameters:
    ----------
    path_to : str or list
        path to input file or files
    fas_name : str, optional
        name of the fasta file
    fas_id : str, optional
        id of the concatenated sequence
    fas_descr : str, optional
        description of the fasta sequence
    """
    if type(path_to) == str:
        cat_seq = _cat_fasta_records(path_to)
    elif type(path_to) == list:
        cat_seq = ""
        for file in path_to:
            cat_seq += _cat_fasta_records(file)
          
    cat_seq.id = fas_id
    cat_seq.description = fas_descr
    SeqIO.write(cat_seq, fas_name, "fasta")        
 

def plot_contigs_cover_gc(path_to):
    """takes spades assembler output which is fasta
    file containing contigs, and
    creates two plots:
    1. distribution of GC content in contigs
    2. GC content vs log2 coverage depth 
    
    Parameters:
    -----------
    path_to : str
        path to input file
    """
    
    container = []
    for seq_record in SeqIO.parse(path_to, "fasta"):
        entry = (seq_record.id, GC(seq_record.seq))
        container.append(entry)
    gc = [x[1] for x in container]
    
    fig = plt.figure()
    sns.distplot(gc, hist=False, kde_kws={"shade":True})
    plt.title("GC_distribution")
    plt.xlabel("GC content, %")
    plt.savefig("contigs_GC_distribution.jpeg", format="jpeg")
    fig.close()
    
    coverage = []
    for el in container:
        cov = el[0].split("_")[-1]
        coverage.append(float(cov))
    cov_log2 = [log2(x) for x in coverage]
    
    fig1 = plt.figure(figsize=(10, 8))
    plt.scatter(gc, cov_log2, s=5)
    plt.xlabel("GC content, %")
    plt.ylabel("log2 coverage depth")
    plt.title("coverage of the contigs vs GC content", fontsize=15)
    plt.savefig("GC_content_vs_contigs_coverage.jpeg", format="jpeg")
    fig1.close()





















      
        
        
        
        
        
        
        
