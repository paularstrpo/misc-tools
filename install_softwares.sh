#!/bin/bash

# ----------------------------------------------------------------------------
#         INSTALL SOFTWARES
#
#         Automate the installation of bioinformatics softwares
#             Paula Restrepo
#             paularstrpo@gmail.com
# ----------------------------------------------------------------------------

ScriptHome=""

mkdir -p ${ScriptHome}/bin
ScriptHome=${ScriptHome}/bin
cd $ScriptHome
# ----------------------------------------------------------------------------
# 1. Burrows-Wheeler Aligner (BWA) - (http://bio-bwa.sourceforge.net/)
# ----------------------------------------------------------------------------

bwa_version="bwa-0.7.17"
curl -L -O https://sourceforge.net/projects/bio-bwa/files/${bwa_version}.tar.bz2

tar -xvf ${bwa_version}.tar.bz2
cd $bwa_version
make

BWA_PATH=${ScriptHome}/$bwa_version
find -name $BWA_PATH -printf

# git clone https://github.com/lh3/bwa.git
# cd bwa; make
# ./bwa index ref.fa
# ./bwa mem ref.fa read-se.fq.gz | gzip -3 > aln-se.sam.gz
# ./bwa mem ref.fa read1.fq read2.fq | gzip -3 > aln-pe.sam.gz

# ----------------------------------------------------------------------------
# 2. GenomeAnalysisTK (GATK) - (https://software.broadinstitute.org/gatk/)
# ----------------------------------------------------------------------------
cd ${ScriptHome}

gatk_version="GenomeAnalysisTK-3.8-0"
curl -L -O https://software.broadinstitute.org/gatk/download/${gatk_version}.tar.bz2

tar -vxjf ${gatk_version}.tar.bz2

mv ${gatk_version}*/GenomeAnalysisTK.jar ${ScriptHome}/${gatk_version}.jar
GATK=${ScriptHome}/${gatk_version}.jar
find -name $GATK -printf

# ----------------------------------------------------------------------------
# 3. Picard Tools - (http://broadinstitute.github.io/picard/)
# ----------------------------------------------------------------------------

cd ${ScriptHome}
picard_version="2.16.0"
curl -L https://github.com/broadinstitute/picard/releases/download/${picard_version}/picard.jar -O ${ScriptHome}/"picard-${picard_version}.jar"
PICARD=${ScriptHome}/"picard-${picard_version}.jar"
find -name $PICARD -printf
