#!/usr/bin/bash

input_file=$1

module load pandoc
module load R

Rscript -e 'library(rmarkdown); rmarkdown::render("'${input_file}'", "html_document")'
