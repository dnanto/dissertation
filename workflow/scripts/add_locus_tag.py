#!/usr/bin/env python

import sys
import re

from Bio import SeqIO

record = SeqIO.read(sys.stdin, "genbank")

key = "locus_tag"
for feature in record.features:
    if feature.type == "CDS" and key not in feature.qualifiers:
        product = feature.qualifiers.get("product", [""])[0]
        if product:
            feature.qualifiers[key] = re.sub(r"[^a-zA-Z0-9*_-]", "_", product)

SeqIO.write(record, sys.stdout, "genbank")
