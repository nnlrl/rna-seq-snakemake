#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Author    : nnlrl
# Created   : 2022/04/10
# Title     : multiqc.py
# Objective :

from os import path

from snakemake.shell import shell


input_dirs = set(path.dirname(fp) for fp in snakemake.input)
output_dir = path.dirname(snakemake.output[0])
output_name = path.basename(snakemake.output[0])
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "multiqc"
    " {snakemake.params}"
    " --force"
    " -o {output_dir}"
    " -n {output_name}"
    " {input_dirs}"
    " {log}"
)
