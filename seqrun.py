# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 11:28:28 2018

@author: babin
"""

import subprocess


###########################################################################
# all the function with _run suffix are the patterns for the wrappers 
# used on the same machine, where you have only to run the function from 
# any script. so, paths to the tools, out and in folders must be
# specified.
    

def bowtie2_run(bowtie_build, ref_seq, prefix, bowtie_align, mapping_name, 
                fastq_file, output, summary):
    """
    runs bowtie2
    
    args: 
    path to bowtie_build.py
    path to reference seq (fasta)
    path to and prefix file name (only id letters)
    name of the run (for log)
    input fastq file
    output .sam file name
    path and name.txt for summary
    """
    
    build_args = [bowtie_build, ref_seq, prefix]
    
    build_inst = subprocess.Popen(build_args, stderr=subprocess.PIPE)
    for line in build_inst.stderr:
        print(line)
        
        
    align_args = [bowtie_align, "-x", prefix, "-U", fastq_file, "-S", output]
    
    align_inst = subprocess.Popen(align_args, stderr=subprocess.PIPE)
    
    with open(summary, "a") as f_obj:
        f_obj.write(mapping_name + "\n")
        for line in align_inst.stderr:
            f_obj.write(str(line) + "\n")
        f_obj.write("***********end of run************")
        
        

def fastqc_run(path_to_fastqc, path_to_file, file_log=False):
    """
    runs fastqc  
    
    'path_to_fastqc' arg should be a file 
    name, path, or a list of names/paths 
    'path_to_file' is path to fastq file
    if file_log=True, stderr of fastqc
    will be written into a file
    """
    
    if type(path_to_file) == list:
        for name_path in path_to_file:
            fastqc_inst = subprocess.Popen([path_to_fastqc, name_path], 
                                           stderr=subprocess.PIPE)
            
            if file_log:
                with open("fastqc_log.txt", "a") as f_obj:
                    for line in fastqc_inst.stderr:
                         f_obj.write(str(line) + "\n")
                    f_obj.write("***************\n")
            
            else:
                for line in fastqc_inst.stderr:
                    print(line)
                
    else:
        fastqc_inst = subprocess.Popen([path_to_fastqc, path_to_file], 
                                       stderr=subprocess.PIPE)
        
        
        # fastqc_inst.kill() # to kill the process
        
        if file_log:
            with open("fastqc_log.txt", "a") as f_obj:
                for line in fastqc_inst.stderr:
                    f_obj.write(str(line) + "\n")
                f_obj.write("***************\n")
        else:
             for line in fastqc_inst.stderr:
                    print(line)



def prefetch_run(sra_ids):
    """
    launch prefetch utility
    
    args:
    path to sra ids in the txt file
    """
    
    with open(sra_ids) as f_obj:
        for line in f_obj:
            line = line.strip()
            prefetch_inst = subprocess.Popen(["/home/yuriy/tools/sratoolkit.2.9.0-centos_linux64/bin/prefetch", line],
                                            stderr=subprocess.PIPE)
            
            for line in prefetch_inst.stderr:
                print(line)
                print("**********************")
            
            


def fastq_dump_run(path_to, out_dir=False):
    """
    run fastq-dump utility
    path to the utility itself is in fastq_dump_inst, 
    change path if needed    

    args:
    path to file with accessions
    path to output folder
    """
    
    if out_dir:
        
        with open(path_to) as f_obj:
            for line in f_obj:
                line = line.strip()
                fastq_dump_inst = subprocess.Popen(["/home/yuriy/tools/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump",
                                                    "--outdir",
                                                    out_dir,
                                                    "--split-files",
                                                    line],
                                                    stderr=subprocess.PIPE)
                for line in fastq_dump_inst.stderr:
                    print(line)
    
    
    else:
        with open(path_to_sra) as f_obj:
            for line in f_obj:
                fastq_dump_inst = subprocess.Popen(["/home/yuriy/tools/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump",
                                                    "--split-files",
                                                    line],
                                                    stderr=subprocess.PIPE)
                for line in fastq_dump_inst.stderr:
                    print(line)
                
         
        
        