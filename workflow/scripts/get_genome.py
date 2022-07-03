#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Author    : nnlrl
# Created   : 2022/01/20
# Title     : get_genome.py
# Objective :

import subprocess
import sys
from itertools import product
from snakemake.shell import shell

species = snakemake.params.species
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

url = f"https://hgdownload.soe.ucsc.edu/goldenPath/{species}/bigZips/{species}.fa.gz"

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
