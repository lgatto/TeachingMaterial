#!/bin/sh

# Test detonate, rsem, and bowtie. Last modified on 2015-02-24 by Brian High

# Works on 64-bit Linux or Mac OSX systems. You must have already installed
# development tools like XCode (OSX). For Ubuntu Linux, you can install with:
# $ sudo apt-get update; sudo apt-get install build-essential
# On a Max (OX) you will also need libpng, which you can install with "brew":
# $ brew install libpng
# You can find out how to install and use brew at: http://brew.sh/
# 
# Run the script as shown below to log standard output and error streams:
# $ time bash -x ~/detonate_test.sh 2>detonate.err 1>detonate.log &
#
# You can view progress with "tail", for example (Ctrl-C to quit):
# $ tail -f detonate.err
# $ tail -f detonate.log

# bowtie must already be installed and found in PATH list. See:
# http://bowtie-bio.sourceforge.net/index.shtml
# For example, to install bowtie for Mac OSX, try this script from our repo:
# $ bash install_bowtie_osx_64.sh   # Read the comments first, though.
# Or on Ubuntu Linux, just run: sudo apt-get install bowtie

# Most of the commands below came from this vignette:
# http://deweylab.biostat.wisc.edu/detonate/vignette.html

# Step 1: Build the prerequisities

# Set PATH to include path to rsem and bowtie binaries.
export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:~/bin

# Check for presence of bowtie in PATH
which bowtie || exit 1

mkdir -p ~/src
cd ~/src

git clone https://github.com/bli25wisc/RSEM.git
cd RSEM/
make
cd ..
export PATH=$PATH:~/src/RSEM

# Set machine type so that blat will compile properly.
# On Linux or Mac OSX systems, the default MACHTYPE  
# (e.g., x86_64-pc-linux-gnu or x86_64-apple-darwin13) 
# will not work with compiling blat. The error you would get is like:
# make[1]: *** No rule to make target `../lib//jkweb.a', 
#     needed by `blat'.  Stop.
export MACHTYPE=$(echo $MACHTYPE | awk -F- '{ print $1 }')

curl -o blatSrc35.zip 'https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip'
unzip blatSrc35.zip 
cd blatSrc/
mkdir -p ~/bin/$MACHTYPE
make
export PATH=$PATH:~/bin/$MACHTYPE
cd ..

# Step 2: Build the DETONATE software

git clone 'https://github.com/deweylab/detonate.git'
cd detonate/
make

# Step 3: Run RSEM-EVAL on each assembly

./rsem-eval/rsem-eval-calculate-score \
    examples/toy_SE.fq \
    examples/toy_assembly_1.fa \
    examples/rsem_eval_1 76 \
    --transcript-length-parameters \
    rsem-eval/true_transcript_length_distribution/mouse.txt -p 16
./rsem-eval/rsem-eval-calculate-score \
    examples/toy_SE.fq \
    examples/toy_assembly_2.fa examples/rsem_eval_2 76 \
    --transcript-length-parameters \
    rsem-eval/true_transcript_length_distribution/mouse.txt -p 16
./rsem-eval/rsem-eval-calculate-score \
    examples/toy_SE.fq \
    examples/toy_assembly_3.fa \
    examples/rsem_eval_3 76 \
    --transcript-length-parameters \
    rsem-eval/true_transcript_length_distribution/mouse.txt -p 16

cat examples/rsem_eval_1.score
cat examples/rsem_eval_1.score | awk '$1 == "Score"'
cat examples/rsem_eval_2.score | awk '$1 == "Score"'
cat examples/rsem_eval_3.score | awk '$1 == "Score"'

# Step 4: Estimate the “true” assembly.

rsem-prepare-reference --bowtie \
    examples/toy_ref.fa examples/toy_rsem_ref
rsem-calculate-expression -p 12 \
    examples/toy_SE.fq examples/toy_rsem_ref examples/toy_rsem_expr

./ref-eval/ref-eval-estimate-true-assembly \
    --reference examples/toy_rsem_ref \
    --expression examples/toy_rsem_expr \
    --assembly examples/ta \
    --alignment-policy best

# Step 5: Compute the kmer-compression score for each assembly.

rsem-prepare-reference --bowtie \
    examples/ta_0.fa examples/ta_0_ref
rsem-calculate-expression -p 12 \
    examples/toy_SE.fq examples/ta_0_ref examples/ta_0_expr

./ref-eval/ref-eval --scores kc \
    --A-seqs examples/toy_assembly_1.fa \
    --B-seqs examples/ta_0.fa \
    --B-expr examples/ta_0_expr.isoforms.results \
    --kmerlen 76 --readlen 76 --num-reads 46988 | tee examples/kc_1.txt
./ref-eval/ref-eval --scores kc \
    --A-seqs examples/toy_assembly_2.fa \
    --B-seqs examples/ta_0.fa \
    --B-expr examples/ta_0_expr.isoforms.results \
    --kmerlen 76 --readlen 76 --num-reads 46988 | tee examples/kc_2.txt
./ref-eval/ref-eval --scores kc \
    --A-seqs examples/toy_assembly_3.fa \
    --B-seqs examples/ta_0.fa \
    --B-expr examples/ta_0_expr.isoforms.results \
    --kmerlen 76 --readlen 76 --num-reads 46988 | tee examples/kc_3.txt

cat examples/kc_1.txt
cat examples/kc_1.txt | awk '$1 == "kmer_compression_score"'
cat examples/kc_2.txt | awk '$1 == "kmer_compression_score"'
cat examples/kc_3.txt | awk '$1 == "kmer_compression_score"'

# Step 6: Compute the alignment-based scores for each assembly.

blat -minIdentity=80 examples/ta_0.fa examples/toy_assembly_1.fa \
    examples/toy_assembly_1_to_ta_0.psl
blat -minIdentity=80 examples/ta_0.fa examples/toy_assembly_2.fa \
    examples/toy_assembly_2_to_ta_0.psl
blat -minIdentity=80 examples/ta_0.fa examples/toy_assembly_3.fa \
    examples/toy_assembly_3_to_ta_0.psl
blat -minIdentity=80 examples/toy_assembly_1.fa examples/ta_0.fa \
   examples/ta_0_to_toy_assembly_1.psl
blat -minIdentity=80 examples/toy_assembly_2.fa examples/ta_0.fa \
   examples/ta_0_to_toy_assembly_2.psl
blat -minIdentity=80 examples/toy_assembly_3.fa examples/ta_0.fa \
    examples/ta_0_to_toy_assembly_3.psl

./ref-eval/ref-eval --scores contig,nucl --weighted no \
    --A-seqs examples/toy_assembly_1.fa \
    --B-seqs examples/ta_0.fa \
    --A-to-B examples/toy_assembly_1_to_ta_0.psl \
    --B-to-A examples/ta_0_to_toy_assembly_1.psl \
    --min-frac-identity 0.90 | tee examples/contig_nucl_1.txt
./ref-eval/ref-eval --scores contig,nucl --weighted no \
    --A-seqs examples/toy_assembly_2.fa \
    --B-seqs examples/ta_0.fa \
    --A-to-B examples/toy_assembly_2_to_ta_0.psl \
    --B-to-A examples/ta_0_to_toy_assembly_2.psl \
    --min-frac-identity 0.90 | tee examples/contig_nucl_2.txt
./ref-eval/ref-eval --scores contig,nucl --weighted no \
    --A-seqs examples/toy_assembly_3.fa \
    --B-seqs examples/ta_0.fa \
    --A-to-B examples/toy_assembly_3_to_ta_0.psl \
    --B-to-A examples/ta_0_to_toy_assembly_3.psl \
    --min-frac-identity 0.90 | tee examples/contig_nucl_3.txt

cat examples/contig_nucl_1.txt
cat examples/contig_nucl_1.txt | awk '$1 == "unweighted_contig_F1"'
cat examples/contig_nucl_2.txt | awk '$1 == "unweighted_contig_F1"'
cat examples/contig_nucl_3.txt | awk '$1 == "unweighted_contig_F1"'
cat examples/contig_nucl_1.txt | awk '$1 == "unweighted_nucl_F1"'
cat examples/contig_nucl_2.txt | awk '$1 == "unweighted_nucl_F1"'
cat examples/contig_nucl_3.txt | awk '$1 == "unweighted_nucl_F1"'
