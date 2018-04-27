#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 12:49:30 2018

@author: yuriy
"""


import subprocess


def prefetch_run(path_to_prefetch, sra_ids):
    """
    launch prefetch utility
    
    args:
    path to prefetch util
    path to sra ids in the txt file
    """
    
    with open(sra_ids) as f_obj:
        for line in f_obj:
            print(line)
            counter = 0
            for ch in line.strip():
                print(counter, ch)
                counter += 1
            line = line.strip()
            for ch in line.strip():
                print(counter, ch)
                counter += 1
            prefetch_inst = subprocess.Popen(["/home/yuriy/tools/sratoolkit.2.9.0-centos_linux64/bin/prefetch", line],
                                            stderr=subprocess.PIPE)
            
            for line in prefetch_inst.stderr:
                print(line)
                print("**********************")
            
            


def fastq_dump_run(path_to_sra, out_dir=False):
    """
    run fastq-dump utility
    
    args:
    path to sra util
    path to output folder
    """
    
    if out_dir:
        
        with open(path_to_sra) as f_obj:
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
                
                
#function call    
# prefetch_run("/home/yuriy/tools/sratoolkit.2.9.0-centos_linux64/bin/prefetch", "/home/yuriy/bos_taurus/SraAccList1.txt")


                   
fastq_dump_run("/home/yuriy/bos_taurus/SraAccList1.txt", out_dir="/home/yuriy/bos_taurus/")