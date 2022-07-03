#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Author    : nnlrl
# Created   : 2022/04/10
# Title     : multiqc.py
# Objective :

from os import path

from snakemake.shell import shell

reads = snakemake.input.reads
GTF_FILE = snakemake.input.gtf_file

outfolder = snakemake.params.get("outfolder", "")
index_dir = snakemake.params.get("index_dir", "")
SALMON_QUANT_EXEC = snakemake.params.get("SALMON_QUANT_EXEC", "")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if(len(snakemake.input.reads) == 1):
    COMMAND = "{SALMON_QUANT_EXEC} -i {index_dir} -l A -p {snakemake.threads} -r {reads} -o {outfolder} --seqBias --gcBias -g {GTF_FILE} {log}"
elif(len(snakemake.input.reads) == 2):
    COMMAND = "{SALMON_QUANT_EXEC} -i {index_dir} -l A -p {snakemake.threads} -1 {reads[0]} -2 {reads[1]} -o {outfolder} --seqBias --gcBias -g {GTF_FILE} {log}"

shell(COMMAND)
