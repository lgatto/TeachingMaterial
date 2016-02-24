#!/bin/sh
# Install bowtie on a Mac OSX 64-bit system. Replace URL with most recent, as needed.
BINDIR="~/bin"
URL='http://hivelocity.dl.sourceforge.net/project/bowtie-bio/bowtie/1.1.1/bowtie-1.1.1-macos-x86_64.zip'
OUTFILE="$(basename "$URL")"
curl -o "$OUTFILE" "$URL"
unzip "$OUTFILE"
cd "$(basename "$OUTFILE" '-macos-x86_64.zip')"
make
mkdir -p "$BINDIR"
cp bowtie* "$BINDIR"/
