# Bioinformatics for Big Omics Data: Using RSEM, a hands-on example
Raphael Gottardo  
February 17, 2015  

## Overview

This short presentation will take you through a short example of RNA quantification using RSEM. You will need to have the following software installed on your machine (Brian should have those setup on the server):

- Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- RSEM (https://github.com/bli25wisc/RSEM)
- sra toolkit for converting sra files to fastq files (http://www.ncbi.nlm.nih.gov/books/NBK158900/)
- Curl for getting fa, gtf, and isoforms files from Dropbox (http://curl.haxx.se/)
- aspera connect for downloading files from the SRA db (optional)
(http://downloads.asperasoft.com/connect2/)
- GNU parallel (optional)
http://www.gnu.org/software/parallel/


## Setup environment

We first prepare the directory structure, these directory are created only if they don't exist, as follows,


```bash
mkdir -p RSEM_test/{GEO,SRA,FASTQ,RSEM,Reference_Genome}
```



```r
# load libraries
library(data.table)
library(reshape2)
library(GEOquery)
library(SRAdb)
```


```r
gd <- getGEO("GSE45735", destdir = "RSEM_test/GEO/")
```

```
## ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45735/matrix/
## Found 1 file(s)
## GSE45735_series_matrix.txt.gz
## Using locally cached version: RSEM_test/GEO//GSE45735_series_matrix.txt.gz
## Using locally cached version of GPL10999 found here:
## RSEM_test/GEO//GPL10999.soft
```

```r
pd <- pData(gd[[1]])
```


```r
if(!file.exists('RSEM_test/SRAmetadb.sqlite')) {
  sqlfile <- getSRAdbFile(destdir = "RSEM_test/")
  }
sra_con <- dbConnect(SQLite(),'RSEM_test/SRAmetadb.sqlite')
sra_tables <- dbListTables(sra_con)
```

## Download SRA files

We now check whether SRA files have been downloaded already, and download the files if they haven't been dowloaded yet


SRA files have been downloaded: TRUE


```r
# Let's just use the first 2 files as an exercise
pd_small <- pd[1:2,]
SRX_number <- rep(NA, nrow(pd_small))
```


```r
# Find an SSH key for Aspera Connect
sshkey <- '/etc/aspera/asperaweb_id_dsa.openssh' # Use this if not found in ~
macsshkey <- '~/Applications/Aspera\\ Connect.app/Contents/Resources/asperaweb_id_dsa.openssh'
if(file.exists(macsshkey)) { sshkey <- macsshkey }
linuxsshkey <- '~/.aspera/connect/etc/asperaweb_id_dsa.openssh'
if(file.exists(linuxsshkey)) { sshkey <- linuxsshkey }

for(i in 1:nrow(pd_small)) {
    gd <- getGEO(pd_small$geo_accession[i], destdir="RSEM_test/GEO/")
    SRX_number[i] <- gsub(".*=SRX", "SRX", gd@header$relation[1])
  
    # Convert to aspera address
	run_accession <- listSRAfile(SRX_number[i], sra_con, fileType = "sra" )$run
	aspera_url <- paste0("anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra", 
                         "/", substr(run_accession,1,3), "/", substr(run_accession,1,6), 
                         "/", run_accession, "/", run_accession, ".sra")
    if(file.exists(sshkey)) {
	    system(paste0('ascp -i ', sshkey ,' -k 1 -T -l200m ', aspera_url, " RSEM_test/SRA"))
    } else { 
        print("Can't find aspera connect sshkey!") 
    }
}
```


```r
## Add the SRX_number information to the pd_small object
pd_small$srx <- SRX_number
```

## Convert SRA files to FASTQ files

We now check whether SRA files have been converted already, and convert them otherwise


SRA files are converted: TRUE


```bash
### Convert to fastq in parallel
parallel -j 2 fastq-dump {} -O RSEM_test/FASTQ/ ::: RSEM_test/SRA/*.sra
```

## Prepare reference genome

We now prepare our reference genome if it hasn't been done yet. 


Reference genome ready: TRUE


```bash
# This is already done, but I leave it here if we need to redo it
cd RSEM_test/Reference_Genome/

# Download reference genome from dropbox.
[ -f 'hg19.fa' ] || curl -L -o 'hg19.fa' 'https://www.dropbox.com/sh/p4cosmsqwtmpce7/AAAP5PQh60sy83GF8i5lFOw0a/hg19.fa?dl=1'
[ -f 'UCSC.gtf' ] || curl -L -o 'UCSC.gtf' 'https://www.dropbox.com/sh/p4cosmsqwtmpce7/AACGDebvGSwK-EFsPa-9f5qEa/UCSC.gtf?dl=1'
[ -f 'knownIsoforms' ] || curl -L -o 'knownIsoforms' 'https://www.dropbox.com/sh/p4cosmsqwtmpce7/AAAfeTr5yKsRmnHd6XW-XhZUa/knownIsoforms?dl=1'

#    Alternatively, you can download the input files from UCSC...
#../../using_rsem_prep_input.sh
#    ...but you will get slightly different results.

rsem-prepare-reference --gtf UCSC.gtf --transcript-to-gene-map knownIsoforms --bowtie2 hg19.fa hg19 
```

## Align FASTQ files to reference genome using RSEM

We now check whether RSEM results exists already, and compute them otherwise

RSEM results are recomputed: TRUE


```bash
# Run in this directory:
cd RSEM_test/RSEM
parallel -j 2 rsem-calculate-expression --bowtie2 -p 4 {} ../Reference_Genome/hg19 {/.} ::: ../FASTQ/*.fastq
```

We can now put all the RSEM counts/TPM results in one matrix, if it hasn't been done yet.


```r
rsem_files <- list.files("RSEM_test/RSEM/", pattern="genes.results", full.names = TRUE)

# Read all files and create a list of data.tables
rsem_list <- lapply(rsem_files, function(x, ...){
    y<-fread(x, drop=c("length", "effective_length", "FPKM")); 
    y$sample_name <- gsub(".genes.results", "", basename(x));
    return(y);})

# Bind all the data.tables to create a long data.table
rsem_data_long <- rbindlist(rsem_list)

# Rename a column, so that it's a bit R friendly
setnames(rsem_data_long, "transcript_id(s)", "transcript_ids")

# Reshape the long data.table to create a matrix
# Assume that missing entries would have a TPM value of 0
rsem_tpm_matrix <- dcast.data.table(rsem_data_long, gene_id+transcript_ids~sample_name, value.var = "TPM", fill=0)
rsem_count_matrix <- dcast.data.table(rsem_data_long, gene_id+transcript_ids~sample_name, value.var = "expected_count", fill=0)

# Keep track of the gene_id - transcript mapping 
rsem_txs_table <- rsem_tpm_matrix[, c("gene_id","transcript_ids"), with=FALSE]
#save(rsem_txs_table, file="RSEM_test/RSEM/rsem_txs_table.Rda")

# Remove the transcript column (we create an annotation table in the next chunk)
rsem_tpm_matrix <- rsem_tpm_matrix[, transcript_ids:=NULL][order(gene_id)]
rsem_count_matrix <- rsem_count_matrix[, transcript_ids:=NULL][order(gene_id)]

# Write to disc
write.csv(rsem_tpm_matrix, file="RSEM_test/RSEM/rsem_tpm_matrix.csv", row.names=FALSE)
write.csv(rsem_count_matrix, file="RSEM_test/RSEM/rsem_count_matrix.csv", row.names=FALSE)
```

### Get sample metadata


```r
# Write the pdata data
write.csv(pd_small, file="RSEM_test/RSEM/rsem_pdata.csv", row.names=FALSE)
```


## Get transcript annotation

Get feature annotations from BioC and create a feature annotation set file.


```r
# Install TxDb packages if you don't already have it.
txdbpkg <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
require(txdbpkg, character.only=TRUE) || {
    source("http://bioconductor.org/biocLite.R")
    biocLite(txdbpkg, suppressUpdates=TRUE, ask=FALSE)
    library(txdbpkg, character.only=TRUE)
}
```

```
## [1] TRUE
```

```r
library(annotate)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# create a table with one transcript per line
#load("RSEM_test/RSEM/rsem_txs_table.Rda")
tx_to_gid <- rsem_txs_table[,list(transcript_id=strsplit(as.character(transcript_ids),",")[[1]]),by="gene_id"]

# map the transcripts to entrez gene ids
tx_to_eid <- data.table(select(txdb, keys = tx_to_gid[,transcript_id], columns="GENEID", keytype="TXNAME"))
setnames(tx_to_eid, c("TXNAME", "GENEID"), c("transcript_id", "entrez_id"))

# Add gene symbol
tx_to_eid[!is.na(entrez_id),gene_symbol:=getSYMBOL(tx_to_eid[!is.na(entrez_id),entrez_id], data='org.Hs.eg')]

# merge all to map information to RSEM gene_ids
tx_table <- merge(tx_to_gid, tx_to_eid, by="transcript_id")
tx_table_reduced <- tx_table[,list(entrez_id=paste(unique(entrez_id),collapse=","), transcript_id=paste(unique(transcript_id),collapse=","), gene_symbol=paste(unique(gene_symbol),collapse=",")),by="gene_id"]

# Write to disc and order by gene_id
write.csv(tx_table_reduced[order(gene_id)], file="RSEM_test/RSEM/rsem_fdata.csv", row.names=FALSE)
```
