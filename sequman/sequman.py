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
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser  # low level fast fasta parser
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import vcf
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


def _format_time_stamp(time_stamp):
    days, day_time = time_stamp.split(" ")
    day_time = day_time.split(":")
    day_time = "_".join(day_time)
    time_stamp = days + "_time-" + day_time
    
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
    # your email in here
    Entrez.email = ""
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
            time_stamp = _format_time_stamp(time_stamp)
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

def _get_fastq_num_records(path_to):
    with open(path_to) as in_handle:
        total_reads = 0
        reads_ids = []
        for title, seq, qual in FastqGeneralIterator(in_handle):
            total_reads += 1
            reads_ids.append(title.split(" ")[0])
        num_uniq_reads = len(set(reads_ids))
        
    return total_reads, num_uniq_reads


def _show_fastq_info(f, total_reads, num_uniq_reads):
    print("file {0} contains:".format(f))
    print("{0} total reads".format(total_reads))
    print("{0} unique reads ids".format(num_uniq_reads))
    print("--------------------------")
    
    

def fastq_info(path_to=False):
    """prints out information about fastq files:
    number of sequences in the file, 
    and number of unique ids in the file

    without arguments takes as an input
    all fastq files in the current dir
    
    Parameters
    ----------
    path_to_fasta : str or list
        path to input file, or list of paths
    """
    
    if type(path_to) == str:
        total_reads, num_uniq_reads = _get_fastq_num_records(path_to)
        _show_fastq_info(path_to, total_reads, num_uniq_reads)
        
    elif type(path_to) == list:
        for path in path_to:
            total_reads, num_uniq_reads = _get_fastq_num_records(path)
            _show_fastq_info(path, total_reads, num_uniq_reads)  
    else:
        current_dir_content = os.listdir()
        for f in current_dir_content:
            if f.rsplit(".", 1)[-1] == "fastq":
                total_reads, num_uniq_reads = _get_fastq_num_records(f)
                _show_fastq_info(f, total_reads, num_uniq_reads)  


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


def cat_fasta_seq(path_to, fas_name="cat_seq.fasta", fas_id="cat_seq", fas_descr=""):
    """concatenates  sequences from fasta files
    into one long sequence. takes one multifasta 
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






def count_indels(vcf_file, min_depth=10, verbose="True"):
    """counts indels in vcf file
    ----------------
    vcf_file: str
        input vcf
    min_depth: int
        minimum depth in favour of indel, 10 by default
    verbose: bool
        True - prints information about the variants
        False - keeps silent
    """

    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    counter = 0

    if verbose:
        for record in vcf_reader:
            if "INDEL" not in record.INFO.keys():
                continue
            elif record.INFO["DP4"][2] + record.INFO["DP4"][3] < min_depth:
                continue
            else:
                print("chromosome: %s, position: %s, ref: %s, indel variant: %s" \
                     % (record.CHROM, record.POS, record.REF, record.ALT ))
                print("depth at position: %s" % record.INFO["DP"])
                print("reads supporting reference: %d" %(record.INFO["DP4"][0] + record.INFO["DP4"][1]))
                print("reads supporting indel variant: %d" %(record.INFO["DP4"][2] + record.INFO["DP4"][3]))
                print("==========================================================================")
                counter += 1
    else:
        for record in vcf_reader:
            if "INDEL" not in record.INFO.keys():
                continue
            elif record.INFO["DP4"][2] + record.INFO["DP4"][3] < min_depth:
                continue
            else:
                counter += 1
    
    print("total number of indels %s" %counter)
    



def count_snps(vcf_file, min_depth=10, verbose="True"):
    """counts SNPs in vcf file
    ----------------
    vcf_file: str
        input vcf
    min_depth: int
        minimum depth in favour of SNPs, 10 by default
    verbose: bool
        True - prints information about the variants
        False - keeps silent
    """

    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    counter = 0

    if verbose:
        for record in vcf_reader:
            if "INDEL" in record.INFO.keys():
                continue
            elif record.INFO["DP4"][2] + record.INFO["DP4"][3] < min_depth:
                continue
            else:
                print("chromosome: %s, position: %s, ref: %s, snp variant: %s" \
                     % (record.CHROM, record.POS, record.REF, record.ALT ))
                print("depth at position: %s" % record.INFO["DP"])
                print("reads supporting reference: %d" %(record.INFO["DP4"][0] + record.INFO["DP4"][1]))
                print("reads supporting snp variant: %d" %(record.INFO["DP4"][2] + record.INFO["DP4"][3]))
                print("==========================================================================")
                counter += 1
    else:
        for record in vcf_reader:
            if "INDEL" in record.INFO.keys():
                continue
            elif record.INFO["DP4"][2] + record.INFO["DP4"][3] < min_depth:
                continue
            else:
                counter += 1
    
    print("total number of SNPs %s" %counter)
    




def vcf_to_df(vcf_file, min_depth=10, var_type="snp"):
    """creates pandas dataframe from the vcf file data
    ----------------
    vcf_file: str
        input vcf
    min_depth: int
        minimum depth in favour of variant, 10 by default
    var_type: str
        snp - prints information about the variants
        indel - keeps silent
    """

    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    
    vcf_data ={"chrom": [], "pos": [], "ref": [], "var": [], "total_depth": [],
              "depth_ref": [], "depth_var": []}
    
    if var_type == "snp":
        for record in vcf_reader:
            if "INDEL" in record.INFO.keys():
                continue
            elif record.INFO["DP4"][2] + record.INFO["DP4"][3] < min_depth:
                continue
            else:
                vcf_data["chrom"].append(record.CHROM)
                vcf_data["pos"].append(record.POS)
                vcf_data["ref"].append(record.REF)
                vcf_data["var"].append(record.ALT)
                vcf_data["total_depth"].append(record.INFO["DP"])
                vcf_data["depth_ref"].append(record.INFO["DP4"][0] + record.INFO["DP4"][1])
                vcf_data["depth_var"].append(record.INFO["DP4"][2] + record.INFO["DP4"][3])
    
    elif var_type == "indel":
        for record in vcf_reader:
            if "INDEL" not in record.INFO.keys():
                continue
            elif record.INFO["DP4"][2] + record.INFO["DP4"][3] < min_depth:
                continue
            else:
                vcf_data["chrom"].append(record.CHROM)
                vcf_data["pos"].append(record.POS)
                vcf_data["ref"].append(record.REF)
                vcf_data["var"].append(record.ALT)
                vcf_data["total_depth"].append(record.INFO["DP"])
                vcf_data["depth_ref"].append(record.INFO["DP4"][0] + record.INFO["DP4"][1])
                vcf_data["depth_var"].append(record.INFO["DP4"][2] + record.INFO["DP4"][3])
    else:
        print("var_type arg not valid")

    df = pd.DataFrame.from_dict(vcf_data)
    
    return df


def _get_ref_seq(input_file):
    """extracts first record from fasta
    which must be a reference
    """
    # low level parser returns tuple of id and sequence
    with open("./" + input_file) as in_handle:
        for title, seq in SimpleFastaParser(in_handle):
            ref_seq = seq
            #ref_seq_id = title
            break
        
    return ref_seq




def _get_input_files():
    valid_extansions = ["fasta", "fas", "fa"]
    input_files = os.listdir("./")
    input_files = [f for f in input_files if f.rsplit(".", 1)[-1] in valid_extansions]
    
    return input_files


def _fill_gaps_from_ref(f, ref_seq, record_container):
    with open(f) as in_handle:
        for title, seq in SimpleFastaParser(in_handle):
            new_seq = ""
            index_counter = 0
            
            for char in seq:
                if char == "-":
                    try:
                        new_seq += ref_seq[index_counter]
                    except IndexError:
                        new_seq += char
                else:
                    new_seq += char
                index_counter += 1
            record_container.append(SeqRecord(Seq(new_seq), id=title,
                                              description=""))
    return record_container


def fill_gaps():
    # to start from command line use: python -c "from sequman import fill_gaps; fill_gaps()"
    input_files = _get_input_files()
    print("job started...")
    for f in input_files:
        ref_seq = _get_ref_seq(f)
        record_container = []
        try:
            record_container = _fill_gaps_from_ref(f, ref_seq, record_container)
        except Exception as e:
            print("error while handling a file: {0}. Error message: {1}".format(f, e))
        
        # here goes writing into a file
        SeqIO.write(record_container, f.rsplit(".",1)[0] + "_gaps_filled.fasta", "fasta")
        print("file {0} done".format(f))
    print("...job done")




def _collect_fasta_records(in_fastas):
    
    record_container = []
    for f in in_fastas:
        with open(f) as in_handle:
            for title, seq in SimpleFastaParser(in_handle):
                record_container.append(SeqRecord(Seq(seq), id=title, description=""))
                
    return record_container



def merge_fasta(in_fastas=False, out_name="appended_fasta_records_"):
    """puts all records from fasta files into a single fasta file
    ----------------
    in_fastas: list or tuple
        input fasta files
    out_name: str
        name for output fasta file
    """
    
    valid_extensions = ["fasta", "fas", "fa"]
    
    if in_fastas:
       record_container = _collect_fasta_records(in_fastas)
        
    else:
        in_fastas = os.listdir("./")
        in_fastas = [f for f in in_fastas if f.rsplit(".", 1)[-1] in valid_extensions]
        record_container = _collect_fasta_records(in_fastas)

    time_stamp = _get_current_time()
    time_stamp = _format_time_stamp(time_stamp)
    
    SeqIO.write(record_container, out_name + time_stamp + ".fasta", "fasta")  
    



def _get_intersected_ids(in_fastas):
          
    ids_container = []
    for f in in_fastas:
        id_set = set()
        with open(f) as in_handle:
            for title, seq in SimpleFastaParser(in_handle):
                id_set.add(title)
            ids_container.append(id_set)
    
    intersected_ids = list(set.intersection(*ids_container))    
    
    return intersected_ids
        


def cat_by_id(in_fastas=False, out_name="cat_by_id_"):
    """
    cats sequences by their ids from separate fasta files 
    into one sequence and outputs a single file
    ----------------
    in_fastas: list or tuple
        input fasta files
    out_name: str
        name for output fasta file
    ----------------    
    NOTE: files must start with number like this : '1_', 
    to denote the order of concatenation
    """
    valid_extensions = ["fasta", "fas", "fa"]
    
    
    if in_fastas:
        intersected_ids = _get_intersected_ids(in_fastas)
    else:
        in_fastas = os.listdir("./")
        in_fastas = [f for f in in_fastas if f.rsplit(".", 1)[-1] in valid_extensions]
        intersected_ids = _get_intersected_ids(in_fastas)
    
    intersected_ids.sort() # to sort record by their number, which goes first
    intersected_records = []
    for seq_id in intersected_ids:
        cat_seq = ""
        for f in in_fastas:
            with open(f) as in_handle:
                for title, seq in SimpleFastaParser(in_handle):
                    if title == seq_id:
                        cat_seq += seq
                        break
        intersected_records.append(SeqRecord(Seq(cat_seq), 
                                                 id=title, 
                                                 description=""))
    
    time_stamp = _get_current_time()
    time_stamp = _format_time_stamp(time_stamp)
    
    SeqIO.write(intersected_records, out_name + time_stamp + ".fasta", "fasta")    



def u_to_t(in_fasta, out_name="changed_u_to_t"):
    """
    replaces U to T in fasta sequences
    saves fasta file
    ----------------
    in_fasta: input fasta files
    out_name: str
        name for output fasta file
    ----------------    
    """
    seq_collection = []
    
    with open(in_fasta) as in_handle:
            for title, seq in SimpleFastaParser(in_handle):
                new_seq = ""
                for char in seq:
                    if char == "U":
                        new_seq += "T"
                    else:
                        new_seq += char
                
                seq_collection.append(SeqRecord(Seq(new_seq), id=title, description=""))
                
                
    SeqIO.write(seq_collection, "out_name", "fasta")  
    
























        
        
        

        
        
        
