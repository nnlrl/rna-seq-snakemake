rule genomeCoverage:
    input:
        bam=os.path.join(MAPPED_READS_DIR, MAPPER, '{sample}_Aligned.sortedByCoord.out.bam'),
        bai=os.path.join(MAPPED_READS_DIR, MAPPER, '{sample}_Aligned.sortedByCoord.out.bam.bai')
    output:
        os.path.join(BIGWIG_DIR, MAPPER, '{sample}.forward.bw'),
        os.path.join(BIGWIG_DIR, MAPPER, '{sample}.reverse.bw'),
        os.path.join(BIGWIG_DIR, MAPPER, '{sample}.bw')
    log:
        os.path.join(LOG_DIR, MAPPER, 'genomeCoverage.forward.{sample}.log'),
        os.path.join(LOG_DIR, MAPPER, 'genomeCoverage.reverse.{sample}.log'),
        os.path.join(LOG_DIR, MAPPER, 'genomeCoverage.{sample}.log')
    resources:
        mem_mb = config["execution"]["rules"]["coverage_bamCoverage"]["memory"]
    conda:
        WORKDIR + "/workflow/envs/bamCoverage.yaml"
    shell:
        """
        {BAMCOVERAGE_EXEC} -b {input.bam} -o {output[0]} --filterRNAstrand forward >> {log[0]} 2>&1
        {BAMCOVERAGE_EXEC} -b {input.bam} -o {output[1]} --filterRNAstrand reverse >> {log[1]} 2>&1
        {BAMCOVERAGE_EXEC} -b {input.bam} --normalizeUsing RPKM --binSize 100 -o {output[2]} >> {log[2]} 2>&1
        """
