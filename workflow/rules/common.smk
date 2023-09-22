import glob

import pandas as pd
from snakemake.remote import FTP
from snakemake.utils import validate

############ vars ############

WORKDIR           = "/home/nnlrl/ngs/rna"

GENOME_FASTA      = config['locations']['genome-fasta']
CDNA_FASTA        = config['locations']['cdna-fasta']
READS_DIR         = config['locations']['reads-dir']
OUTPUT_DIR        = config['locations']['output-dir']
ORGANISM          = config['organism']
MAPPER            = config['mapping']['mapper']
GENOME_BUILD      = config['mapping']['genome_build']

SCRIPTS_DIR       = WORKDIR + '/workflow/scripts'

TRIMMED_READS_DIR = os.path.join(OUTPUT_DIR, 'trimmed_reads')
LOG_DIR           = os.path.join(OUTPUT_DIR, 'logs')
QC_DIR            = os.path.join(OUTPUT_DIR, 'QC')
MULTIQC_DIR       = os.path.join(OUTPUT_DIR, 'multiqc')
MAPPED_READS_DIR  = os.path.join(OUTPUT_DIR, 'mapped_reads')
BIGWIG_DIR        = os.path.join(OUTPUT_DIR, 'bigwig_files')
COUNTS_DIR        = os.path.join(OUTPUT_DIR, 'feature_counts')
SALMON_DIR        = os.path.join(OUTPUT_DIR, 'salmon_output')


def toolArgs(name):
    if 'args' in config['tools'][name]:
        return config['tools'][name]['args']
    else:
        return ""

def tool(name):
    cmd = config['tools'][name]['executable']
    return cmd + " " + toolArgs(name)


MULTIQC_EXEC       = tool('multiqc')
STAR_EXEC_MAP      = tool('star_map')
STAR_EXEC_INDEX    = tool('star_index')
HISAT2_EXEC        = tool('hisat2')
HISAT2_BUILD_EXEC  = tool('hisat2-build')
SALMON_INDEX_EXEC  = tool('salmon_index')
SALMON_QUANT_EXEC  = tool('salmon_quant')
SAMTOOLS_EXEC      = tool('samtools')
GUNZIP_EXEC        = tool('gunzip') # for STAR
RSCRIPT_EXEC       = tool('Rscript')
SED_EXEC           = tool('sed')
FASTP_EXEC         = tool('fastp')
BAMCOVERAGE_EXEC   = tool('bamCoverage')

STAR_INDEX_THREADS    = config['execution']['rules']['star_index']['threads']
HISAT2_BUILD_THREADS  = config['execution']['rules']['hisat2-build']['threads']
HISAT2_THREADS        = config['execution']['rules']['hisat2']['threads']
STAR_MAP_THREADS      = config['execution']['rules']['star_map']['threads']
SALMON_INDEX_THREADS  = config['execution']['rules']['salmon_index']['threads']
SALMON_QUANT_THREADS  = config['execution']['rules']['salmon_quant']['threads']
FEATURECOUNTS_THREADS = config['execution']['rules']['featureCounts']['threads']

GTF_FILE          = config['locations']['gtf-file']
SAMPLE_SHEET_FILE = config['locations']['sample-sheet']

DE_ANALYSIS_LIST  = config.get('DEanalyses', {})


for analysis in DE_ANALYSIS_LIST.keys():
    DE_ANALYSIS_LIST[analysis]['covariates'] = (
        DE_ANALYSIS_LIST[analysis]['covariates'] if 'covariates' in DE_ANALYSIS_LIST[analysis].keys()
        else ''
    )

## Load sample sheet
with open(SAMPLE_SHEET_FILE, 'r') as fp:
  rows =  [row for row in csv.reader(fp, delimiter=',')]
  header = rows[0]; rows = rows[1:]
  SAMPLE_SHEET = [dict(zip(header, row)) for row in rows]

# Convenience function to access fields of sample sheet columns that
# match the predicate.  The predicate may be a string.
def lookup(column, predicate, fields=[]):
    if inspect.isfunction(predicate):
        records = [line for line in SAMPLE_SHEET if predicate(line[column])]
    else:
        records = [line for line in SAMPLE_SHEET if line[column]==predicate]
    return [record[field] for record in records for field in fields]

SAMPLES = [line['name'] for line in SAMPLE_SHEET]

############ vars ############




############ targets ############

targets = {
    # rule to print all rule descriptions
    'help': {
        'description': "Print all rules and their descriptions.",
        'files': []
    },
    'final-report': {
        'description': "Produce a comprehensive report.  This is the default target.",
        'files':
        [os.path.join(MULTIQC_DIR, 'multiqc_report.html'),
         os.path.join(COUNTS_DIR, "raw_counts", "salmon", "counts_from_SALMON.transcripts.tsv"),
         os.path.join(COUNTS_DIR, "raw_counts", "salmon", "counts_from_SALMON.genes.tsv"),
         os.path.join(COUNTS_DIR, "normalized", "salmon", "TPM_counts_from_SALMON.transcripts.tsv"),
         os.path.join(COUNTS_DIR, "normalized", "salmon", "TPM_counts_from_SALMON.genes.tsv"),
         os.path.join(COUNTS_DIR, "raw_counts", MAPPER, "counts.tsv"),
         os.path.join(COUNTS_DIR, "normalized", MAPPER, "deseq_normalized_counts.tsv"),
         os.path.join(COUNTS_DIR, "normalized", MAPPER, "deseq_size_factors.txt")] +
        expand(os.path.join(BIGWIG_DIR, MAPPER, '{sample}.forward.bw'), sample = SAMPLES) +
        expand(os.path.join(BIGWIG_DIR, MAPPER, '{sample}.reverse.bw'), sample = SAMPLES) +
        expand(os.path.join(BIGWIG_DIR, MAPPER, '{sample}.bw'), sample = SAMPLES) +
        expand(os.path.join(OUTPUT_DIR, "report", MAPPER, '{analysis}.deseq.report.html'), analysis = DE_ANALYSIS_LIST.keys()) +
        expand(os.path.join(OUTPUT_DIR, "report", 'salmon', '{analysis}.salmon.transcripts.deseq.report.html'), analysis = DE_ANALYSIS_LIST.keys()) +
        expand(os.path.join(OUTPUT_DIR, "report",  'salmon', '{analysis}.salmon.genes.deseq.report.html'), analysis = DE_ANALYSIS_LIST.keys())
    },
    'deseq_report_star': {
        'description': "Produce one HTML report for each analysis based on STAR results.",
        'files':
          expand(os.path.join(OUTPUT_DIR, "report", 'star', '{analysis}.deseq.report.html'), analysis = DE_ANALYSIS_LIST.keys())
    },
    'deseq_report_hisat2': {
        'description': "Produce one HTML report for each analysis based on Hisat2 results.",
        'files':
          expand(os.path.join(OUTPUT_DIR, "report", 'hisat2', '{analysis}.deseq.report.html'), analysis = DE_ANALYSIS_LIST.keys())
    },
    'deseq_report_salmon_transcripts': {
        'description': "Produce one HTML report for each analysis based on SALMON results at transcript level.",
        'files':
          expand(os.path.join(OUTPUT_DIR, "report", 'salmon', '{analysis}.salmon.transcripts.deseq.report.html'), analysis = DE_ANALYSIS_LIST.keys())
    },
    'deseq_report_salmon_genes': {
        'description': "Produce one HTML report for each analysis based on SALMON results at gene level.",
        'files':
          expand(os.path.join(OUTPUT_DIR, "report", "salmon", '{analysis}.salmon.genes.deseq.report.html'), analysis = DE_ANALYSIS_LIST.keys())
    },
    'star_map' : {
        'description': "Produce a STAR mapping results in BAM file format.",
        'files':
          expand(os.path.join(MAPPED_READS_DIR, "star", '{sample}_Aligned.sortedByCoord.out.bam'), sample = SAMPLES)
    },
    'star_counts': {
        'description': "Get count matrix from STAR mapping results using summarizeOverlaps.",
        'files':
          [os.path.join(COUNTS_DIR, "raw_counts", "star", "counts.tsv")]
    },
    'hisat2_map' : {
        'description': "Produce Hisat2 mapping results in BAM file format.",
        'files':
          expand(os.path.join(MAPPED_READS_DIR, "hisat2", '{sample}_Aligned.sortedByCoord.out.bam'), sample = SAMPLES)
    },
    'hisat2_counts': {
        'description': "Get count matrix from Hisat2 mapping results using summarizeOverlaps.",
        'files':
          [os.path.join(COUNTS_DIR, "raw_counts", "hisat2", "counts.tsv")]
    },
    'genome_coverage': {
        'description': "Compute genome coverage values from BAM files - save in bigwig format",
        'files':
          expand(os.path.join(BIGWIG_DIR, MAPPER, '{sample}.forward.bw'), sample = SAMPLES) +
          expand(os.path.join(BIGWIG_DIR, MAPPER, '{sample}.reverse.bw'), sample = SAMPLES) +
          expand(os.path.join(BIGWIG_DIR, MAPPER, '{sample}.bw'), sample = SAMPLES)
    },
    'salmon_index' : {
        'description': "Create SALMON index file.",
        'files':
          [os.path.join(OUTPUT_DIR, 'salmon_index', "pos.bin")]
    },
    'salmon_quant' : {
        'description': "Calculate read counts per transcript using SALMON.",
        'files':
          expand(os.path.join(SALMON_DIR, "{sample}", "quant.sf"), sample = SAMPLES) +
	  expand(os.path.join(SALMON_DIR, "{sample}", "quant.genes.sf"), sample = SAMPLES)
    },
    'salmon_counts': {
        'description': "Get count matrix from SALMON quant.",
        'files':
          [os.path.join(COUNTS_DIR, "raw_counts", "salmon", "counts_from_SALMON.transcripts.tsv"),
	   os.path.join(COUNTS_DIR, "raw_counts", "salmon", "counts_from_SALMON.genes.tsv"),
	   os.path.join(COUNTS_DIR, "normalized", "salmon", "TPM_counts_from_SALMON.transcripts.tsv"),
	   os.path.join(COUNTS_DIR, "normalized", "salmon", "TPM_counts_from_SALMON.genes.tsv")]
    },
    'multiqc': {
        'description': "Get multiQC report based on alignments and QC reports.",
        'files':
          [os.path.join(MULTIQC_DIR, 'multiqc_report.html')]
    }
}

############ targets ############


# Selected output files from the above set.
selected_targets = config['execution']['target'] or ['final-report']

# FIXME: the list of files must be flattened twice(!).  We should make
# sure that the targets really just return simple lists.
from itertools import chain
OUTPUT_FILES = list(chain.from_iterable([targets[name]['files'] for name in selected_targets]))


# determine if the sample library is single end or paired end
def isSingleEnd(args):
    sample = args[0]
    files = lookup('name', sample, ['reads', 'reads2'])
    count = sum(1 for f in files if f)
    if count == 2:
        return False
    elif count == 1:
        return True


# function to pass read files to trim/filter/qc improvement
def trim_reads_input(args):
    sample = args[0]
    return [os.path.join(READS_DIR, f) for f in lookup('name', sample, ['reads', 'reads2']) if f]


def map_input(args):
    sample = args[0]
    reads_files = [os.path.join(READS_DIR, f) for f in lookup('name', sample, ['reads', 'reads2']) if f]
    if len(reads_files) > 1:
        return [os.path.join(TRIMMED_READS_DIR, "{sample}.trimmed.R1.fq.gz".format(sample=sample)), os.path.join(TRIMMED_READS_DIR, "{sample}.trimmed.R2.fq.gz".format(sample=sample))]
    elif len(reads_files) == 1:
        return [os.path.join(TRIMMED_READS_DIR, "{sample}.trimmed.fq.gz".format(sample=sample))]


# I cannot do function composition, so it's gotta be this awkward definition instead.
def hisat2_file_arguments(args):
    files = map_input(args)
    if len(files) == 2:
        return "-1 {} -2 {}".format(files[0], files[1])
    elif len(files) == 1:
        return "-U {}".format(files[0])
