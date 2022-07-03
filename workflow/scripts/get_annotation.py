#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Author    : nnlrl
# Created   : 2022/01/20
# Title     : get_annotation.py
# Objective :

import subprocess
import sys
from snakemake.shell import shell

# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz

species = snakemake.params.species
fmt = snakemake.params.fmt
assembly = snakemake.params.assembly


log = snakemake.log_fmt_shell(stdout=False, stderr=True)

suffix = ""
if fmt == "gtf":
    suffix = "gtf.gz"
elif fmt == "gff3":
    suffix = "gff3.gz"

url = f"https://hgdownload.soe.ucsc.edu/goldenPath/{species}/bigZips/genes/{species}.{assembly}.{suffix}"

try:
    shell("(curl --retry 100 --retry-delay 1 -L --url {url} | gzip -d > {snakemake.output[0]}) {log}")
except subprocess.CalledProcessError as e:
    if snakemake.log:
        sys.stderr = open(snakemake.log[0], "a")
    print(
        "Unable to download annotation data from UCSC. "
        "Did you check that this combination of species, build, and release is actually provided?",
        file=sys.stderr,
    )
    exit(1)
