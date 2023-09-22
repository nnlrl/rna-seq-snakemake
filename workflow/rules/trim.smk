# fastp both trims/filters reads and outputs QC reports in html/json format
rule trim_qc_reads_pe:
    input:
        trim_reads_input
    output:
        r1=os.path.join(TRIMMED_READS_DIR, "{sample}.trimmed.R1.fq.gz"),
        r2=os.path.join(TRIMMED_READS_DIR, "{sample}.trimmed.R2.fq.gz"),
        html=os.path.join(QC_DIR, "{sample}.pe.fastp.html"),
        json=os.path.join(QC_DIR, "{sample}.pe.fastp.json") #notice that multiqc recognizes files ending with fast.json
    log:
        os.path.join(LOG_DIR, 'trim_reads.{sample}.log')
    conda:
        WORKDIR + "/workflow/envs/fastp.yaml"
    shell:
        "{FASTP_EXEC} --in1 {input[0]} --in2 {input[1]} --out1 {output.r1} --out2 {output.r2} -h {output.html} -j {output.json} >> {log} 2>&1"


# fastp both trims/filters reads and outputs QC reports in html/json format
rule trim_qc_reads_se:
    input:
        trim_reads_input
    output:
        r = os.path.join(TRIMMED_READS_DIR, "{sample}.trimmed.fq.gz"),
        html=os.path.join(QC_DIR, "{sample}.se.fastp.html"),
        json=os.path.join(QC_DIR, "{sample}.se.fastp.json") #notice that multiqc recognizes files ending with fast.json
    log:
        os.path.join(LOG_DIR, 'trim_reads.{sample}.log')
    conda:
        WORKDIR + "/workflow/envs/fastp.yaml"
    shell:
        "{FASTP_EXEC} --in1 {input[0]} --out1 {output.r} -h {output.html} -j {output.json} >> {log} 2>&1 "
