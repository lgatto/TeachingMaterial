#!/bin/bash

# This script prepares the three input files needed for Using_RSEM.Rmd.
# Or you can download them (hg19.fa, UCSC.gtf, and knownIsoforms) from:
# https://www.dropbox.com/sh/p4cosmsqwtmpce7/AADgXoeK4RwBTfDlDUr_dwQYa?dl=0
# We do this in the "Prepare reference genome" section of Using_RSEM.Rmd.
#
# Parts of this script were adapted from:
#     http://watson.nci.nih.gov/~sdavis/tutorials/biowulf-2011/
#     http://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format

# Required utilities: 
# - bash which basename curl tar cat touch grep chmod gunzip genePredToGtf
#
# Notes:
# - genePredToGtf is from kentUtils: https://github.com/ENCODE-DCC/kentUtils
# - The rest will be standard on a Unix, Linux, or Mac OSX system.
# - This script has not been tested on a Windows/Cygwin system. 

# Steps performed by this script are:
# 1. Download and combine the human genome hg19 chromosome "fa" files
# 2. Download known isoforms file
# 3. Download known gene file and convert to GTF file

# After successfully running this script, you should be able to run:
#    rsem-prepare-reference --gtf UCSC.gtf \
#        --transcript-to-gene-map knownIsoforms \ 
#        --bowtie2 hg19.fa hg19

# Make sure we have the tools we need
for i in curl tar touch grep chmod cat gunzip genePredToGtf; do \
    which "$i" >/dev/null
    if [ $? -ne 0 ]; then \
        echo "Can't find $i. Aborting!" && exit 1
    fi
done

# Set the parent URL hosting the human genome "hg19" data files
PARENTURL='http://hgdownload.cse.ucsc.edu/goldenPath/hg19'

# Download, extract and combine chromosome files (unless already done)
URL="$PARENTURL/bigZips/chromFa.tar.gz"
BIGFA="hg19.fa"
[ -s "$BIGFA" ] || (curl -s -o - "$URL" | tar xzf - -O > "$BIGFA")

# Download and extract known isoforms file (unless already done)
URL="$PARENTURL/database/knownIsoforms.txt.gz"
ISOFORMS=$(basename "$URL" .txt.gz)
[ -s "$ISOFORMS" ] || (curl -s "$URL" | gunzip -c - > "$ISOFORMS")

# Download known gene file and create GTF file
GTF='UCSC.gtf'
HGCONF=~/'.hg.conf'
touch "$HGCONF"
chmod 600 "$HGCONF"
grep -q 'db.host=genome-mysql.cse.ucsc.edu' "$HGCONF" || echo '
db.host=genome-mysql.cse.ucsc.edu
db.user=genomep
db.password=password
central.db=hgcentral' >> "$HGCONF"
[ -s "$GTF" ] || genePredToGtf hg19 knownGene "$GTF"
