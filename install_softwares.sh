#!/bin/bash

# ----------------------------------------------------------------------------
#         INSTALL SOFTWARES
#
#         Automate the installation of bioinformatics softwares
#             Paula Restrepo
#             paularstrpo@gmail.com
# ----------------------------------------------------------------------------


mkdir -p $HOME/bin
ScriptHome=$HOME/bin
cd $ScriptHome

# ----------------------------------------------------------------------------
# Versions
# ----------------------------------------------------------------------------


bwa_version="bwa-0.7.17"
samtools_version="samtools-1.6"
bcftools_version="bcftools-1.6"
htslib_version="htslib-1.6"
gatk_version="GenomeAnalysisTK-3.8-0"
picard_version="2.16.0"


# ----------------------------------------------------------------------------
# GATK & PICARD
# ----------------------------------------------------------------------------

curl -L -O https://software.broadinstitute.org/gatk/download/${gatk_version}.tar.bz2
curl -L https://github.com/broadinstitute/picard/releases/download/${picard_version}/picard.jar -O ${ScriptHome}/"picard-${picard_version}.jar"
tar -vxjf ${gatk_version}.tar.bz2
mv ${gatk_version}*/GenomeAnalysisTK.jar ${ScriptHome}/${gatk_version}.jar

exit

# ----------------------------------------------------------------------------
# Download packages
# ----------------------------------------------------------------------------


curl -L -O https://sourceforge.net/projects/bio-bwa/files/${bwa_version}.tar.bz2
curl -L -O https://github.com/samtools/samtools/releases/download/1.6/${samtools_version}.tar.bz2
curl -L -O https://github.com/samtools/bcftools/releases/download/1.6/${bcftools_version}.tar.bz2
curl -L -O https://github.com/samtools/htslib/releases/download/1.6/${htslib_version}.tar.bz2

tar -xvf ${bwa_version}.tar.bz2
tar -vxjf ${samtools_version}.tar.bz2
tar -vxjf ${bcftools_version}.tar.bz2
tar -vxjf ${htslib_version}.tar.bz2

rm *.tar.*

# ----------------------------------------------------------------------------
# BWA
# ----------------------------------------------------------------------------


cd $bwa_version
make

#find -name $BWA_PATH -type d -printf

# git clone https://github.com/lh3/bwa.git
# cd bwa; make
# ./bwa index ref.fa
# ./bwa mem ref.fa read-se.fq.gz | gzip -3 > aln-se.sam.gz
# ./bwa mem ref.fa read1.fq read2.fq | gzip -3 > aln-pe.sam.gz

# ----------------------------------------------------------------------------
# Samtools, BCFtools & HTSlib - (http://www.htslib.org/)
# ----------------------------------------------------------------------------

./configure --prefix=${samtools_version}
make

./configure --prefix=${bcftools_version}
make

./configure --prefix=${htslib_version}
make



# ----------------------------------------------------------------------------
# Homebrew (MacOS only)
# ----------------------------------------------------------------------------
#/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

# ----------------------------------------------------------------------------
# Circos (MacOS/local machine only)
# ----------------------------------------------------------------------------
#
# #Install requisite perl libraries
# brew install cpanminus
# cpanm Config::General Font::TTF::Font Math::Bezier Math::VecStat Readonly Set::IntSpan Text::Format Statistics::Basic SVG
#
# mkdir ${ScriptHome}/lib ; cd ${ScriptHome}/lib
# curl -O http://wangqinhu.com/data/gd/gd.tar.gz
# tar -zxf gd.tar.gz
# cd gd
# sudo ./install
#
# #Install Circos
# mkdir ${ScriptHome}/circos ; cd ${ScriptHome}/circos
# circos_pkg=''
# circos_tools=''
# circos_tut=''
#
# curl -O http://circos.ca/distribution/${circos_pkg}.tgz
# curl -O http://circos.ca/distribution/${circos_tools}.tgz
# curl -O http://circos.ca/distribution/${circos_tut}.tgz
#
# tar xvfz $circos_pkg
# tar xvfz $circos_tools
# tar xvfz $circos_tut
#
# ln -s $circos_pkg current
