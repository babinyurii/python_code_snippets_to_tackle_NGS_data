#### a bunch of functions for routine work with sequence data:

`fetch_seq()` downloads sequences from nucleotide database by id nums and saves them in genbank format

`blast_fasta()` blast records from a fasta file, writes results into a tab-delimited txt file

`fasta_info()` prints out information about fasta files: number of sequences in the file, sequence id numbers, lengths of sequences and GC content

`fastq_info()` prints out information about fastq files: number of sequences in the file, and number of unique ids in the file

`split_fasta()` splits fasta file containing several sequences into the corresponding number of fasta files

`cat_fasta_seq()` oncatenates  sequences from fasta files into one long sequence. takes one multifasta  or several fasta files as an input

`plot_contigs_cover_gc()`takes spades assembler output which is fasta file containing contigs, and creates two plots:  distribution of GC content in contigs, and GC content vs log2 coverage depth 

`count_indels()` counts indels in a vcf file

`count_snps()` counts SNPs in a vcf file

`vcf_to_df()` creates pandas dataframe from the vcf file data

`merge_fasta()` puts all records from input fasta files into a single fasta file

`cat_by_id()` concatenates sequences by their id from separate fasta files into a single sequence. outputs fasta with concatenated sequences
