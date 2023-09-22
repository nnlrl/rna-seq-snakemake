rule counts_from_SALMON:
    input:
        quantFiles = expand(os.path.join(SALMON_DIR, "{sample}", "quant.sf"), sample=SAMPLES),
        quantGenesFiles = expand(os.path.join(SALMON_DIR, "{sample}", "quant.genes.sf"), sample=SAMPLES),
        colDataFile = rules.translate_sample_sheet_for_report.output
    output:
        os.path.join(COUNTS_DIR, "raw_counts", "salmon", "counts_from_SALMON.transcripts.tsv"),
        os.path.join(COUNTS_DIR, "raw_counts", "salmon","counts_from_SALMON.genes.tsv"),
        os.path.join(COUNTS_DIR, "normalized", "salmon", "TPM_counts_from_SALMON.transcripts.tsv"),
        os.path.join(COUNTS_DIR, "normalized", "salmon", "TPM_counts_from_SALMON.genes.tsv")
    resources:
        mem_mb = config["execution"]["rules"]["count_reads"]["memory"]
    log:
        os.path.join(LOG_DIR, "salmon", 'salmon_import_counts.log')
    conda:
        WORKDIR + "/workflow/envs/Rscript.yaml"
    shell:
        "{RSCRIPT_EXEC} {SCRIPTS_DIR}/counts_matrix_from_SALMON.R {SALMON_DIR} {COUNTS_DIR} {input.colDataFile} >> {log} 2>&1"


# rule count_reads:
#     input:
#         bam = os.path.join(MAPPED_READS_DIR, MAPPER, "{sample}_Aligned.sortedByCoord.out.bam"),
#         bai = os.path.join(MAPPED_READS_DIR, MAPPER, "{sample}_Aligned.sortedByCoord.out.bam.bai")
#     output:
#         os.path.join(MAPPED_READS_DIR, MAPPER, "{sample}.read_counts.csv")
#     resources:
#         mem_mb = 5000
#     log:
#         os.path.join(LOG_DIR, MAPPER, "{sample}.count_reads.log")
#     params:
#         single_end   = isSingleEnd,
#         mode         = config['counting']['counting_mode'],
#         nonunique    = config['counting']['drop_nonunique'],
#         strandedness = config['counting']['strandedness'],
#         feature      = config['counting']['feature'],
#         group_by     = config['counting']['group_feature_by'],
#         yield_size   = config['counting']['yield_size']
#     conda:
#         "../envs/Rscript.yaml"
#     shell:
#         "{RSCRIPT_EXEC} {SCRIPTS_DIR}/count_reads.R {wildcards.sample} {input.bam} {GTF_FILE} {params.single_end} {params.mode} {params.nonunique} {params.strandedness} {params.feature} {params.group_by} {params.yield_size} >> {log} 2>&1"


rule featureCounts:
    input:
        samples = os.path.join(MAPPED_READS_DIR, MAPPER, "{sample}_Aligned.sortedByCoord.out.bam"),
        annotation = GTF_FILE
    output:
        os.path.join(MAPPED_READS_DIR, MAPPER, "{sample}.read_counts.csv")
    resources:
        mem_mb = config["execution"]["rules"]["featureCounts"]["memory"]
    log:
        os.path.join(LOG_DIR, "{sample}.featureCounts.log")
    params:
        single_end   = isSingleEnd,
        strandedness = config['counting']['strandedness'],
        feature      = config['counting']['feature'],
        group_by     = config['counting']['group_feature_by'],
        annotation_file_type = config['counting']['annotation_file_type'],
    conda:
        WORKDIR + "/workflow/envs/featureCounts.yaml"
    threads:
        FEATURECOUNTS_THREADS
    script:
        WORKDIR + "/workflow/scripts/featureCounts.py"


rule collate_read_counts:
    input:
        expand(os.path.join(MAPPED_READS_DIR, MAPPER, "{sample}.read_counts.csv"), sample = SAMPLES)
    output:
        os.path.join(COUNTS_DIR, "raw_counts", MAPPER, "counts.tsv")
    resources:
        mem_mb = 1024
    log:
        os.path.join(LOG_DIR, MAPPER, "collate_read_counts.log")
    params:
        mapped_dir = os.path.join(MAPPED_READS_DIR, MAPPER),
        script = os.path.join(SCRIPTS_DIR, "collate_read_counts.R")
    conda:
        WORKDIR + "/workflow/envs/Rscript.yaml"
    shell:
        "{RSCRIPT_EXEC} {params.script} {params.mapped_dir} {output} >> {log} 2>&1"


# create a normalized counts table including all samples
# using the median-of-ratios normalization procedure of
# deseq2
rule norm_counts_deseq:
    input:
        counts_file = os.path.join(COUNTS_DIR, "raw_counts", MAPPER, "counts.tsv"),
        colDataFile = rules.translate_sample_sheet_for_report.output
    output:
        size_factors = os.path.join(COUNTS_DIR, "normalized", MAPPER, "deseq_size_factors.txt"),
        norm_counts = os.path.join(COUNTS_DIR, "normalized", MAPPER, "deseq_normalized_counts.tsv")
    resources:
        mem_mb = 1024
    log:
        os.path.join(LOG_DIR, MAPPER, "norm_counts_deseq.log")
    params:
        script=os.path.join(SCRIPTS_DIR, "norm_counts_deseq.R"),
        outdir=os.path.join(COUNTS_DIR, "normalized", MAPPER)
    conda:
        WORKDIR + "/workflow/envs/Rscript.yaml"
    shell:
        "{RSCRIPT_EXEC} {params.script} {input.counts_file} {input.colDataFile} {params.outdir} >> {log} 2>&1"


rule report1:
    input:
        counts=os.path.join(COUNTS_DIR, "raw_counts", MAPPER, "counts.tsv"),
        coldata=str(rules.translate_sample_sheet_for_report.output),
    params:
        outdir=os.path.join(OUTPUT_DIR, "report", MAPPER),
        reportR=os.path.join(SCRIPTS_DIR, "runDeseqReport.R"),
        reportRmd=os.path.join(SCRIPTS_DIR, "deseqReport.Rmd"),
        case = lambda wildcards: DE_ANALYSIS_LIST[wildcards.analysis]['case_sample_groups'],
        control = lambda wildcards: DE_ANALYSIS_LIST[wildcards.analysis]['control_sample_groups'],
        covariates = lambda wildcards: DE_ANALYSIS_LIST[wildcards.analysis]['covariates'],
    log:
        os.path.join(LOG_DIR, MAPPER, "{analysis}.report.log")
    output:
        os.path.join(OUTPUT_DIR, "report", MAPPER, '{analysis}.deseq.report.html')
    resources:
        mem_mb = config["execution"]["rules"]["reports"]["memory"]
    conda:
        WORKDIR + "/workflow/envs/Rscript.yaml"
    shell:
        "{RSCRIPT_EXEC} {params.reportR} --prefix='{wildcards.analysis}' --reportFile={params.reportRmd} --countDataFile={input.counts} --colDataFile={input.coldata} --gtfFile={GTF_FILE} --caseSampleGroups='{params.case}' --controlSampleGroups='{params.control}' --covariates='{params.covariates}'  --workdir={params.outdir} --organism='{ORGANISM}'  >> {log} 2>&1"

rule report2:
    input:
        counts=os.path.join(COUNTS_DIR, "raw_counts", "salmon", "counts_from_SALMON.transcripts.tsv"),
        coldata=str(rules.translate_sample_sheet_for_report.output)
    params:
        outdir=os.path.join(OUTPUT_DIR, "report", 'salmon'),
        reportR=os.path.join(SCRIPTS_DIR, "runDeseqReport.R"),
        reportRmd=os.path.join(SCRIPTS_DIR, "deseqReport.Rmd"),
        case = lambda wildcards: DE_ANALYSIS_LIST[wildcards.analysis]['case_sample_groups'],
        control = lambda wildcards: DE_ANALYSIS_LIST[wildcards.analysis]['control_sample_groups'],
        covariates = lambda wildcards: DE_ANALYSIS_LIST[wildcards.analysis]['covariates'],
    log:
        os.path.join(LOG_DIR, "salmon", "{analysis}.report.salmon.transcripts.log")
    output:
        os.path.join(OUTPUT_DIR, "report", 'salmon', '{analysis}.salmon.transcripts.deseq.report.html')
    resources:
        mem_mb = config["execution"]["rules"]["reports"]["memory"]
    conda:
        WORKDIR + "/workflow/envs/Rscript.yaml"
    shell:
        "{RSCRIPT_EXEC} {params.reportR} --prefix='{wildcards.analysis}.salmon.transcripts' --reportFile={params.reportRmd} --countDataFile={input.counts} --colDataFile={input.coldata} --gtfFile={GTF_FILE} --caseSampleGroups='{params.case}' --controlSampleGroups='{params.control}' --covariates='{params.covariates}' --workdir={params.outdir} --organism='{ORGANISM}' >> {log} 2>&1"

rule report3:
    input:
        counts=os.path.join(COUNTS_DIR, "raw_counts", "salmon", "counts_from_SALMON.genes.tsv"),
        coldata=str(rules.translate_sample_sheet_for_report.output)
    params:
        outdir=os.path.join(OUTPUT_DIR, "report", "salmon"),
        reportR=os.path.join(SCRIPTS_DIR, "runDeseqReport.R"),
        reportRmd=os.path.join(SCRIPTS_DIR, "deseqReport.Rmd"),
        case = lambda wildcards: DE_ANALYSIS_LIST[wildcards.analysis]['case_sample_groups'],
        control = lambda wildcards: DE_ANALYSIS_LIST[wildcards.analysis]['control_sample_groups'],
        covariates = lambda wildcards: DE_ANALYSIS_LIST[wildcards.analysis]['covariates'],
    log:
        os.path.join(LOG_DIR, "salmon", "{analysis}.report.salmon.genes.log")
    output:
        os.path.join(OUTPUT_DIR, "report", "salmon", '{analysis}.salmon.genes.deseq.report.html')
    resources:
        mem_mb = config["execution"]["rules"]["reports"]["memory"]
    conda:
        WORKDIR + "/workflow/envs/Rscript.yaml"
    shell:
        "{RSCRIPT_EXEC} {params.reportR} --prefix='{wildcards.analysis}.salmon.genes' --reportFile={params.reportRmd} --countDataFile={input.counts} --colDataFile={input.coldata} --gtfFile={GTF_FILE} --caseSampleGroups='{params.case}' --controlSampleGroups='{params.control}' --covariates='{params.covariates}' --workdir={params.outdir} --organism='{ORGANISM}' >> {log} 2>&1"
