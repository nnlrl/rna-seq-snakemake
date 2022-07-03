#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Author    : nnlrl
# Created   : 2022/04/10
# Title     : samtools-faidx.py
# Objective :

from snakemake.shell import shell
from snakemake_wrapper_utils.samtools import get_samtools_opts

samtools_opts = get_samtools_opts(
    snakemake, parse_threads=False, parse_write_index=False, parse_output_format=False
)
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell("samtools faidx {samtools_opts} {extra} {snakemake.input[0]} {log}")
