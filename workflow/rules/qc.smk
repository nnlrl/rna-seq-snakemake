rule multiqc:
    input:
        salmon_output=expand(os.path.join(SALMON_DIR, "{sample}", "quant.sf"), sample = SAMPLES),
        mapping_output=expand(os.path.join(MAPPED_READS_DIR, MAPPER, '{sample}_Aligned.sortedByCoord.out.bam'), sample=SAMPLES)
    output:
        os.path.join(MULTIQC_DIR, 'multiqc_report.html')
    resources:
        mem_mb = 200
    log:
        os.path.join(LOG_DIR, f'multiqc.{MAPPER}.log')
    conda:
        WORKDIR + "/workflow/envs/multiqc.yaml"
    shell:
        "{MULTIQC_EXEC} -f -o {MULTIQC_DIR} {OUTPUT_DIR} >> {log} 2>&1"
