rule star_index:
    input:
        GENOME_FASTA
    output:
        star_index_file = os.path.join(OUTPUT_DIR, 'star_index', "SAindex")
    resources:
        mem_mb = config["execution"]["rules"]["star_index"]["memory"]
    params:
        star_index_dir = os.path.join(OUTPUT_DIR, 'star_index')
    log:
        os.path.join(LOG_DIR, 'star_index.log')
    conda:
        WORKDIR + "/workflow/envs/star.yaml"
    threads:
        STAR_INDEX_THREADS
    shell:
        "{STAR_EXEC_INDEX} --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.star_index_dir} --genomeFastaFiles {input} --sjdbGTFfile {GTF_FILE} >> {log} 2>&1"


rule hisat2_index:
    input:
        GENOME_FASTA
    output:
        [os.path.join(OUTPUT_DIR, "hisat2_index", f"{GENOME_BUILD}_index.{n}.ht2l") for n in [1, 2, 3, 4, 5, 6, 7, 8]]
    resources:
        mem_mb = config["execution"]["rules"]["hisat2-build"]["memory"]
    params:
        index_directory = os.path.join(OUTPUT_DIR, "hisat2_index"),
    log:
        os.path.join(LOG_DIR, 'hisat2_index.log')
    conda:
        WORKDIR + "/workflow/envs/hisat2.yaml"
    threads:
        HISAT2_BUILD_THREADS
    shell:
        "{HISAT2_BUILD_EXEC} -p {threads} --large-index {input} {params.index_directory}/{GENOME_BUILD}_index >> {log} 2>&1"


rule salmon_index:
    input:
        CDNA_FASTA
    output:
        salmon_index_file = os.path.join(OUTPUT_DIR, 'salmon_index', "hash.bin")
    resources:
        mem_mb = config["execution"]["rules"]["salmon_index"]["memory"]
    params:
        salmon_index_dir = os.path.join(OUTPUT_DIR, 'salmon_index')
    log:
        os.path.join(LOG_DIR, "salmon", 'salmon_index.log')
    conda:
        WORKDIR + "/workflow/envs/salmon.yaml"
    threads:
        SALMON_INDEX_THREADS
    shell:
        "{SALMON_INDEX_EXEC} -t {input} -i {params.salmon_index_dir} -p {threads} >> {log} 2>&1"



rule star_map:
    input:
        # This rule really depends on the whole directory (see
        # params.index_dir), but we can't register it as an input/output
        # in its own right since Snakemake 5.
        index_file = rules.star_index.output.star_index_file,
        reads = map_input
    output:
        os.path.join(MAPPED_READS_DIR, 'star', '{sample}_Aligned.sortedByCoord.out.bam')
    resources:
        mem_mb = config["execution"]["rules"]["star_map"]["memory"]
    params:
        index_dir = rules.star_index.params.star_index_dir,
        output_prefix=os.path.join(MAPPED_READS_DIR, 'star', '{sample}_')
    log:
        os.path.join(LOG_DIR, 'star', 'star_map_{sample}.log')
    conda:
        WORKDIR + "/workflow/envs/star.yaml"
    threads:
        STAR_MAP_THREADS
    shell:
        "{STAR_EXEC_MAP} --runThreadN {threads} --genomeDir {params.index_dir} --readFilesIn {input.reads} --readFilesCommand '{GUNZIP_EXEC} -c' --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.output_prefix} >> {log} 2>&1"


rule hisat2_map:
    input:
        index_files = rules.hisat2_index.output,
        reads = map_input
    output:
        os.path.join(MAPPED_READS_DIR, 'hisat2', '{sample}_Aligned.sortedByCoord.out.bam')
    resources:
        mem_mb = config["execution"]["rules"]["hisat2"]["memory"]
    params:
        samfile = lambda wildcards: os.path.join(MAPPED_READS_DIR, 'hisat2', "_".join([wildcards.sample, 'Aligned.out.sam'])),
        index_dir = rules.hisat2_index.params.index_directory,
        args = hisat2_file_arguments
    log:
        os.path.join(LOG_DIR, 'hisat2', 'hisat2_map_{sample}.log'),
        os.path.join(LOG_DIR, 'hisat2', 'samtools.hisat2.{sample}.log')
    conda:
        WORKDIR + "/workflow/envs/hisat2.yaml"
    threads:
        HISAT2_THREADS
    shell:
        """
        {HISAT2_EXEC} -x {params.index_dir}/{GENOME_BUILD}_index -p {threads} -q -S {params.samfile} {params.args} >> {log[0]} 2>&1
        {SAMTOOLS_EXEC} view -bh {params.samfile} | {SAMTOOLS_EXEC} sort -o {output} >> {log[1]} 2>&1
        rm {params.samfile}
        """


rule index_bam:
    input:
        os.path.join(MAPPED_READS_DIR, MAPPER, '{sample}_Aligned.sortedByCoord.out.bam')
    output:
        os.path.join(MAPPED_READS_DIR, MAPPER, '{sample}_Aligned.sortedByCoord.out.bam.bai')
    resources:
        mem_mb = config["execution"]["rules"]["index_bam"]["memory"]
    conda:
        WORKDIR + "/workflow/envs/samtools.yaml"
    shell:
        "{SAMTOOLS_EXEC} index {input} {output}"


rule salmon_quant:
    input:
        # This rule really depends on the whole directory (see
        # params.index_dir), but we can't register it as an input/output
        # in its own right since Snakemake 5.
        index_file = rules.salmon_index.output.salmon_index_file,
        reads = map_input,
        gtf_file = GTF_FILE
    output:
        os.path.join(SALMON_DIR, "{sample}", "quant.sf"),
        os.path.join(SALMON_DIR, "{sample}", "quant.genes.sf")
    resources:
        mem_mb = config["execution"]["rules"]["salmon_quant"]["memory"]
    params:
        index_dir = rules.salmon_index.params.salmon_index_dir,
        outfolder = os.path.join(SALMON_DIR, "{sample}"),
        SALMON_QUANT_EXEC = {SALMON_QUANT_EXEC},
    log:
        os.path.join(LOG_DIR, "salmon", 'salmon_quant_{sample}.log')
    conda:
        WORKDIR + "/workflow/envs/salmon.yaml"
    threads:
        SALMON_QUANT_THREADS
    script:
        WORKDIR + "/workflow/scripts/salmon_quant.py"
