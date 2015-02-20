#!/bin/bash

# This script prepares the three input files needed for Using _RSEM.Rmd.
# Or you can download them (hg19.fa, UCSC.gtf, and knownIsoforms) from:
# https://www.dropbox.com/sh/p4cosmsqwtmpce7/AADgXoeK4RwBTfDlDUr_dwQYa?dl=0
# We do this in the "Prepare reference genome" section of Using_RSEM.Rmd.

# Required utilities: 
#     bash curl tar touch grep chmod ls mv rm gunzip genePredToGtf (kentUtils)

# Parts of this script were adapted from:
#     http://watson.nci.nih.gov/~sdavis/tutorials/biowulf-2011/
#     http://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format

# Steps performed by this script are:
# 1. Download the human genome hg19 chromosome "fa" files (chromFa.tar.gz)
# 2. Combine the individual "fa" files into a single "hg19.fa" file
# 3. Download known isoforms file
# 4. Create GTF file from known gene file

# Make sure we have the tools we need
for i in curl tar touch grep chmod ls mv rm echo gunzip genePredToGtf; do \
    which $i >/dev/null
    if [ $? -ne 0 ]; then \
        echo "Can't find $i. Aborting!" && exit 1
    fi
done

# Initialize variables
HG19PATH='http://hgdownload.cse.ucsc.edu/goldenPath/hg19'
CHROMEPATH='bigZips/chromFa.tar.gz'
URL="$HG19PATH/$CHROMEPATH"
TARBALL=$(basename "$URL")
FA="hg19.fa"

# Download the chromosome files and extract (unless already done)
[ -s "$TARBALL" ] || curl -O "$URL"
ls chr*.fa 2>/dev/null  1>/dev/null || tar xvzf "$TARBALL"

# Combine the main chromosome files into one genome file
[ -s "$FA" ] || (for i in $(seq 22) X Y M; do \
    cat "chr${i}.fa" >> "$FA"
done)
rm -f chr*.fa

# Download and extract known isoforms file (unless already done)
ISOFORMSPATH='database/knownIsoforms.txt.gz'
URL="$HG19PATH/$ISOFORMSPATH"
TARBALL=$(basename "$URL")
ISOFORMS=$(basename "$TARBALL" .gz)
[ -s "$ISOFORMS" ] || (curl -O "$URL" && gunzip "$TARBALL")
[ -s "$ISOFORMS" ] && mv "$ISOFORMS" $(basename "$ISOFORMS" .txt)

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
#
# Alternatively, you can follow the steps from the URL below to get the GTF:
# https://groups.google.com/a/soe.ucsc.edu/d/msg/genome/kyk7AAm4R-M/9LkE-CRjzioJ
# However, you should note the following quote from the RSEM README: 
#     https://github.com/bli25wisc/RSEM#-usage 
# "Please note that GTF files generated from the UCSC Table Browser do not 
# contain isoform-gene relationship information. However, if you use the UCSC 
# Genes annotation track, this information can be recovered by downloading the 
# knownIsoforms.txt file for the appropriate genome."

# Now you should be ready to run rsem-prepare-reference on these files:
#     rsem-prepare-reference --gtf UCSC.gtf 
#         --transcript-to-gene-map knownIsoforms 
#         --bowtie2 hg19.fa hg19
