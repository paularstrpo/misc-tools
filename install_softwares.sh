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
#find -name $BWA_PATH -type d -printf

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
#find -name $GATK -printf

# ----------------------------------------------------------------------------
# 3. Picard Tools - (http://broadinstitute.github.io/picard/)
# ----------------------------------------------------------------------------
cd ${ScriptHome}
picard_version="2.16.0"
curl -L https://github.com/broadinstitute/picard/releases/download/${picard_version}/picard.jar -O ${ScriptHome}/"picard-${picard_version}.jar"
tar -vxjf ${picard_version}.tar.bz2
PICARD=${ScriptHome}/"picard-${picard_version}.jar"
#find -name $PICARD -printf

# ----------------------------------------------------------------------------
# 4. Samtools, BCFtools & HTSlib - (http://www.htslib.org/)
# ----------------------------------------------------------------------------
cd ${ScriptHome}
samtools_version="samtools-1.6"
bcftools_version="bcftools-1.6"
htslib_version="htslib-1.6"


curl -L -O https://github.com/samtools/samtools/releases/download/1.6/${samtools_version}.tar.bz2
curl -L -O https://github.com/samtools/bcftools/releases/download/1.6/${bcftools_version}.tar.bz2
curl -L -O https://github.com/samtools/htslib/releases/download/1.6/${htslib_version}.tar.bz2

tar -vxjf ${samtools_version}.tar.bz2
tar -vxjf ${bcftools_version}.tar.bz2
tar -vxjf ${htslib_version}.tar.bz2


./configure --prefix=${samtools_version}
make

./configure --prefix=${bcftools_version}
make

./configure --prefix=${htslib_version}
make


# ----------------------------------------------------------------------------
# 5. Homebrew (MacOS only)
# ----------------------------------------------------------------------------
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

# ----------------------------------------------------------------------------
# 6. Circos (MacOS version)
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
