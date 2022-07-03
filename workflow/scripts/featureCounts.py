#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# Author    : nnlrl
# Created   : 2022/06/19
# Title     : featureCounts.py
# Objective :

import tempfile
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
# extra = snakemake.params.get("extra", "")

# optional input files and directories
strand = snakemake.params.get("strandedness", 0)
if int(strand) not in [0, 1, 2]:
    print("Acceptable strandedness options are '0(unspecific), 1(forward), or 2(reverse).")
    exit(1)

annotation_file_type = snakemake.params.get("annotation_file_type", "")
if annotation_file_type:
    annotation_file_type = f"-F {annotation_file_type}"

group_feature_by = snakemake.params.get("group_by", "")
if group_feature_by:
    group_feature_by = f"-g {group_feature_by}"

feature = snakemake.params.get("feature", "")
if feature:
    feature = f"-t {feature}"


singleEnd = snakemake.params.get("single_end", False)
if singleEnd:
    isPair = ""
else:
    isPair = "-p"


with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "featureCounts"
        " -T {snakemake.threads}"
        " -s {strand}"
        " -a {snakemake.input.annotation}"
        " {isPair}"
        " {annotation_file_type}"
        " {group_feature_by}"
        " {feature}"
        " --tmpDir {tmpdir}"
        " -o {snakemake.output[0]}"
        " {snakemake.input.samples}"
        " {log}"
    )
