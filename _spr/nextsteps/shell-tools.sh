grep '^>' EcoliK12_enolase_UPSsimga_NB.fasta | cut -d'|' -f1 | sed 's/>//g' | sort | uniq | grep -v tr| grep -v sp | wc
