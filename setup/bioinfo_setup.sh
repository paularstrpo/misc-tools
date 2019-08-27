#!/usr/bin/bash

conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda update -n base conda
conda create -n py27 python=2.7 anaconda
conda install bamtools bcftools blast blat bowtie2 circexplorer2 beagle circos cnvkit cufflinks entrez-direct emboss fastq-tools fastqc freebayes gatk gatk4 ggplot hisat2 htseq htslib igv igvtools java-jdk kraken multiqc pysam tabix uniprot varscan vcfanno bedtools vcftools bowtie star bwa