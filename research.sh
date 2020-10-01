#!/usr/bin/env bash

# database

blastdb="blast"
mkdir -p "$blastdb"

awk 'NR > 1 { print $1 }' genome.tsv | \
  xargs -L 250 | tr ' ' ',' | \
  xargs -L 1 efetch -db nuccore -format fasta -id > "$blastdb/genome.fasta"

# prototype

awk 'NR > 1 && $0 !~ /^#/ && $6 != "gen" {print $0}' query.tsv | \
while IFS=$'\t' read -r id seq_start seq_stop strand species region protein_id
do
  root="${species}-${region}"
  mkdir -p "$root"
  efetch -db nuccore -id "$id" -seq_start "$seq_start" -seq_stop "$seq_stop" -strand "$strand" -format gb > "$root/ref.gb"
  efetch -db nuccore -id "$id" -seq_start "$seq_start" -seq_stop "$seq_stop" -strand "$strand" -format fasta > "$root/ref.fasta"
  glsearch36 -T 16 -m 8C "$root/ref.fasta" "$blastdb/genome.fasta" > "$root/hits.tsv"
done
