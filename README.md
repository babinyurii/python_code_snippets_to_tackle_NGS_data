`count_snp()` counts SNPs in a vcf file with given coverage depth

`count_indels()` counts indels in a vcf file with given coverage depth

`vcf_to_df()` creates pandas dataframe using vcf file data

`fetch_seq()` downloads sequences in fasta or genbank format from NCBI nucleotide database

`fasta_info()` prints out information about fasta sequences: number of records within a file, their length and GC content

`split_fasta()` splits fasta file containing multiple records into the corresponding number of fasta files

`cat_fasta()` concatenates sequences from one multifasta file or several files into one sequence and writes it as fasta

`plot_contigs_cover_gc()` takes SPAdes assembler output (fasta file containing contigs) and creates two plots: GC distribution and GC content vs contigs coverage:

![gc_vs_coverage](plot_examples/contigs_GC_distribution.jpeg)

![gc_vs_coverage](plot_examples/GC_content_vs_contigs_coverage.jpeg)


`blast_fasta()` blasts sequences from fasta file using nucleotide blas and saves results as a tab-delimited text file

`fastq_info()` prints out information about fastq files: total number of reads in a file and number of unique reads ids
