#!/usr/bin/env bash

# database

blastdb="blast"
mkdir -p "$blastdb"

awk 'NR > 1 { print $1 }' genome.tsv > "$blastdb/genome.ssv"

awk 'NR > 1 { print $1 }' genome.tsv | \
xargs -L 250 | \
tr ' ' ',' | \
xargs -L 1 efetch -db nuccore -format fasta -id | \
makeblastdb \
  -input_type fasta -dbtype nucl -title genome \
  -parse_seqids -hash_index -out "$blastdb/genome" -blastdb_version 5 \
  -logfile "$blastdb/genome.log" -taxid_map "$blastdb/genome.ssv"

# prototype

awk 'NR > 1 && $0 !~ /^#/ {print $0}' query.tsv | \
while IFS=$'\t' read -r id seq_start seq_stop strand species region
do
  root="${species}-${region}"
  mkdir -p "$root"
  efetch -db nuccore -id "$id" -seq_start "$seq_start" -seq_stop "$seq_stop" -strand "$strand" -format gb > "$root/ref.gb"
  efetch -db nuccore -id "$id" -seq_start "$seq_start" -seq_stop "$seq_stop" -strand "$strand" -format fasta >  "$root/ref.fasta"
  blastn \
    -task blastn \
    -db "$blastdb/genome" \
    -query "$root/ref.fasta" \
    -outfmt "7 std qlen sstrand" \
    -num_alignments 100000 \
    -num_threads 16 \
    -out "$root/hits.tsv"
done
