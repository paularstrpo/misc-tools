#!/bin/sh

#This Script is designed to install Circos for WSL (Windows Subsystem for Linux).
#It may also work on debian-based linux systems, but this has not been tested as of 08.08.2016
#Note that you will have to modify your PATH independently of this script, preferably using ~/.bashrc or some equivalent

#Insert the path of where you want to install Circos-- 
	#NOTES:
	#1. For WSL, navigate there using Bash and use PWD
	#
	#2. For WSL, also make sure to use a path under /mnt/c/ which is where your Windows file system should be-- 
	#   otherwise you won't be able to access your circos files through the Windows file explorer!!

PPath='/mnt/c/PATHTOFOLDER' # paste HERE
echo 'your installation path will be: ' ${PPath}'software/circos'

#paste the links to the circos program, tools, and tutorials tarballs
#in each respective variable (between the quotes)

srcLINK = 'http://circos.ca/distribution/circos-0.69-3.tgz'
toolsLINK = 'http://circos.ca/distribution/circos-tools-0.22.tgz'
tutorialsLINK = 'http://circos.ca/distribution/circos-tutorials-0.67.tgz'


#names for output program folders from download-- modifyable by user (optional) 
# *** MUST HAVE A .tgz EXTENSION***
srcNAME = 'circos-0.69-3.tgz'
toolsNAME = 'circos-tools.tgz'
tutorialsNAME = 'circos-tutorials.tgz'


##	##	##

#Make the folder and install circos
cd $PPath
pwd
mkdir $PPath/software
mkdir $PPath/software/circos
cd $PPath/software/circos
pwd

echo 'fetching Circos program and files...'
wget -O $srcNAME $srcLINK
wget -O $toolsNAME $toolsLINK
wget -O $tutorialsNAME $tutorialsLINK
echo 'now to unpack the tarballs...'
tar xvfz $tutorialsNAME
tar xvfz $toolsNAME
tar xvfz $srcNAME
rm *.tgz

$srcNAME = ${srcNAME%/.tgz/}
$toolsNAME = ${toolsNAME%/.tgz/}
$tutorialsNAME = ${tutorialsNAME%/.tgz/}

ln -s $srcNAME current
echo 'download & extraction DONE'
cd $PPath

#Install the required perl modules:
	#NOTE: I had issues installing perl modules using CPAN in WSL
	#      but apt-get seems to work fine so I'll use that here
echo -e \
"installing required perl modules...
-----------------------------------------------------
WARNING: this installation may not be complete: 
after installing, make sure to check for missing modules using circos -modules
-----------------------------------------------------

starting....

"

sudo apt-get -y install \
		libconfig-general-perl \
		libreadonly-perl \
		libfont-ttf-perl \
		libgd-perl \
		liblist-moreutils-perl \
		libmath-bezier-perl \
		libmath-round-perl \
		libmath-vecstat-perl \
		libparams-validate-perl \
		libregexp-common-perl \
		libset-intspan-perl \
		libtext-format-perl \
		libtext-formattable-perl \
		libtext-wrapper-perl \
		libconfig-tiny-perl \
		libconfig-merge-perl \
		libgd-svg-perl \
		libsvg-perl \
		liborlite-statistics-perl \
		libstatistics-basic-perl \
		libclone-perl \
		libcgi-application-plugin-config-simple-perl

sudo apt-get update
sudo apt-get -y upgrade

echo -e \
"-----------------------------------------------------
if Config::General doesn't work then try going into circos-XX.XX/lib/Modules.pm and remove the version number from Config::General
-----------------------------------------------------

-----------------------------------------------------
WARNING: make sure to append the following line (as is) to your ~/.bashrc or equivalent:

export PATH=$PPath/software/circos/current/bin:\$PATH
-----------------------------------------------------

DONE"