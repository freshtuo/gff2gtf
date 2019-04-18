# gff2gtf
convert gff annotation file to gtf format

The purpose of the script is to convert NCBI virus gff format gene annotation file to gtf format, for running STAR, HTSeq-counts, etc

As a result, only 'exon' entries and related entries (parental of 'exon') are collected and reformated into a gtf file, other entries (e.g. 'intron') are discarded. The output gtf file including the needed information ('gene_id', 'transcript_id', etc).

