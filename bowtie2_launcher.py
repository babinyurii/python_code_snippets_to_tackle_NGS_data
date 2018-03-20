#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import subprocess

def bowtie2_run(bowtie_build, ref_seq, prefix, bowtie_align, mapping_name, 
                fastq_file, output, summary):
    """runs bowtie2
    
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
        

bowtie2_run("/home/yuriy/tools/bowtie2-2.3.4.1-linux-x86_64/bowtie2-build",
            "/home/yuriy/sal/minor_fraction/CP008928.1 Salmonella enterica subsp. enterica serovar Enteritidis strain SEJ.fasta",
            "/home/yuriy/sal/salmonella_SEJ",
            "/home/yuriy/tools/bowtie2-2.3.4.1-linux-x86_64/bowtie2",
            "test_mapping",
            "/home/yuriy/sal/sal_data/sal_raw_004.fastq",
            "/home/yuriy/sal/bowtie_out/salmonella_SEJ.sam",
            "/home/yuriy/sal/bowtie_test_log.txt")           
            
    
    
    
    
    
