#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Author    : nnlrl
# Created   : 2022/04/10
# Title     : cutadapt-pe.py
# Objective :

from snakemake.shell import shell


n = len(snakemake.input)
assert n == 2, "Input must contain 2 (paired-end) elements."

extra = snakemake.params.get("extra", "")
adapters = snakemake.params.get("adapters", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

assert (
    extra != "" or adapters != ""
), "No options provided to cutadapt. Please use 'params: adapters=' or 'params: extra='."

shell(
    "cutadapt"
    " {adapters}"
    " {extra}"
    " -o {snakemake.output.fastq1}"
    " -p {snakemake.output.fastq2}"
    " -j {snakemake.threads}"
    " {snakemake.input}"
    " > {snakemake.output.qc} {log}"
)
