#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 12:30:29 2018
@author: yuriy
"""

import subprocess


def fastqc_run(path_tofastqc, name_path_tofile, file_log=False):
    """runs fastqc  
    
    'name_path_tofile' arg should be a file 
    name, path, or a list of names/paths 
    if file_log=True, stderr of fastqc
    will be written into a file
    """
    
    if type(name_path_tofile) == list:
        for name_path in name_path_tofile:
            fastqc_inst = subprocess.Popen([path_tofastqc, name_path], 
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
        fastqc_inst = subprocess.Popen([path_tofastqc, name_path_tofile], 
                                       stderr=subprocess.PIPE)
        
        if file_log:
            with open("fastqc_log.txt", "a") as f_obj:
                for line in fastqc_inst.stderr:
                    f_obj.write(str(line) + "\n")
                f_obj.write("***************\n")
        else:
             for line in fastqc_inst.stderr:
                    print(line)
                
            
# function call    
fastqc_run("/usr/local/bin/fastqc", 
              ["/home/yuriy/sal/sal_raw_003.fastq",
               "/home/yuriy/sal/sal_raw_004.fastq",
               "/home/yuriy/sal/sal_raw_014.fastq",
               "/home/yuriy/sal/sal_raw_015.fastq"], 
               file_log=True)
        

