#!/usr/bin/env bash


# prototype

efetch -db nuccore -id NC_001460.1 -format fasta > A-gen.fasta
efetch -db nuccore -id NC_001460.1 -format genbank > A-gen.gb
efetch -db nuccore -id NC_001460.1 -seq_start 4953 -seq_stop 8513 -format fasta -revcomp > A-pol.fasta
efetch -db nuccore -id NC_001460.1 -seq_start 13394 -seq_stop 14887 -format fasta > A-pen.fasta
efetch -db nuccore -id NC_001460.1 -seq_start 17740 -seq_stop 20499 -format fasta > A-hex.fasta
efetch -db nuccore -id NC_001460.1 -seq_start 29368 -seq_stop 31131 -format fasta > A-fib.fasta

efetch -db nuccore -id NC_011203.1 -format fasta > B-gen.fasta
efetch -db nuccore -id NC_011203.1 -format gb > B-gen.gb
efetch -db nuccore -id NC_011203.1 -seq_start 5051 -seq_stop 8623 -format fasta -revcomp > B-pol.fasta
efetch -db nuccore -id NC_011203.1 -seq_start 13905 -seq_stop 15539 -format fasta > B-pen.fasta
efetch -db nuccore -id NC_011203.1 -seq_start 18418 -seq_stop 21252 -format fasta > B-hex.fasta
efetch -db nuccore -id NC_011203.1 -seq_start 31368 -seq_stop 32327 -format fasta > B-fib.fasta

efetch -db nuccore -id NC_001405.1 -format fasta > C-gen.fasta
efetch -db nuccore -id NC_001405.1 -format gb > C-gen.gb
efetch -db nuccore -id NC_001405.1 -seq_start 5187 -seq_stop 8774 -format fasta -revcomp > C-pol.fasta
efetch -db nuccore -id NC_001405.1 -seq_start 14151 -seq_stop 15866 -format fasta > C-pen.fasta
efetch -db nuccore -id NC_001405.1 -seq_start 18838 -seq_stop 21744 -format fasta > C-hex.fasta
efetch -db nuccore -id NC_001405.1 -seq_start 31030 -seq_stop 32778 -format fasta > C-fib.fasta

efetch -db nuccore -id NC_010956.1 -format fasta > D-gen.fasta
efetch -db nuccore -id NC_010956.1 -format gb > D-gen.gb
efetch -db nuccore -id NC_010956.1 -seq_start 5002 -seq_stop 8523 -format fasta -revcomp > D-pol.fasta
efetch -db nuccore -id NC_010956.1 -seq_start 13516 -seq_stop 15075 -format fasta > D-pen.fasta
efetch -db nuccore -id NC_010956.1 -seq_start 17790 -seq_stop 20621 -format fasta > D-hex.fasta
efetch -db nuccore -id NC_010956.1 -seq_start 30914 -seq_stop 32002 -format fasta > D-fib.fasta

efetch -db nuccore -id NC_003266.2 -format fasta > E-gen.fasta
efetch -db nuccore -id NC_003266.2 -format gb > E-gen.gb
efetch -db nuccore -id NC_003266.2 -seq_start 5033 -seq_stop 8605 -format fasta -revcomp > E-pol.fasta
efetch -db nuccore -id NC_003266.2 -seq_start 13815 -seq_stop 15422 -format fasta > E-pen.fasta
efetch -db nuccore -id NC_003266.2 -seq_start 18248 -seq_stop 21058 -format fasta > E-hex.fasta
efetch -db nuccore -id NC_003266.2 -seq_start 31649 -seq_stop 32926 -format fasta > E-fib.fasta

efetch -db nuccore -id NC_001454.1 -format fasta > F-gen.fasta
efetch -db nuccore -id NC_001454.1 -format gb > F-gen.gb
efetch -db nuccore -id NC_001454.1 -seq_start 4750 -seq_stop 8307 -format fasta -revcomp > F-pol.fasta
efetch -db nuccore -id NC_001454.1 -seq_start 13241 -seq_stop 14755 -format fasta > F-pen.fasta
efetch -db nuccore -id NC_001454.1 -seq_start 17643 -seq_stop 20414 -format fasta > F-hex.fasta
efetch -db nuccore -id NC_001454.1 -seq_start 28751 -seq_stop 29914 -format fasta > F-fib.fasta

efetch -db nuccore -id NC_006879.1 -format fasta > G-gen.fasta
efetch -db nuccore -id NC_006879.1 -format gb > G-gen.gb
efetch -db nuccore -id NC_006879.1 -seq_start 4702 -seq_stop 8247 -format fasta -revcomp > G-pol.fasta
efetch -db nuccore -id NC_006879.1 -seq_start 13129 -seq_stop 14640 -format fasta > G-pen.fasta
efetch -db nuccore -id NC_006879.1 -seq_start 17513 -seq_stop 20308 -format fasta > G-hex.fasta
efetch -db nuccore -id NC_006879.1 -seq_start 28731 -seq_stop 29822 -format fasta > G-fib1.fasta
efetch -db nuccore -id NC_006879.1 -seq_start 29855 -seq_stop 31537 -format fasta > G-fib2.fasta

# genome

awk 'NR > 1 { print $1 }' genome.tsv > genome.ssv

awk 'NR > 1 { print $1 }' genome.tsv | \
xargs -L 250 | \
tr ' ' ',' | \
xargs -L 1 efetch -db nuccore -format fasta -id | \
makeblastdb \
  -input_type fasta -dbtype nucl -title genome \
  -parse_seqids -hash_index -out genome -blastdb_version 5 \
  -logfile genome.log -taxid_map genome.ssv

cat [A-G]-*.fasta | \
blastn \
  -task blastn \
  -db genome \
  -outfmt "7 std qlen sstrand" \
  -num_alignments 10000 \
  -num_threads 16 \
  -out hits.tsv
