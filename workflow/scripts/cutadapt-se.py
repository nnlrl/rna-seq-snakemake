#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Author    : nnlrl
# Created   : 2022/04/10
# Title     : cutadapt-se.py
# Objective :

from snakemake.shell import shell


n = len(snakemake.input)
assert n == 1, "Input must contain 1 (single-end) element."

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
    " -j {snakemake.threads}"
    " -o {snakemake.output.fastq}"
    " {snakemake.input[0]}"
    " > {snakemake.output.qc} {log}"
)
