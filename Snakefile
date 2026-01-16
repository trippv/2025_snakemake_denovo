import pandas as pd

# Load configuration file
configfile: "config/config.yaml"

# Read samples table
SAMPLES_TABLE = pd.read_csv(config["samples_file"], sep="\t")

# Filter included samples
INCLUDED_SAMPLES = SAMPLES_TABLE[SAMPLES_TABLE["include"] == 1]

# Dictionary with sample data
SAMPLE_DICT = {
    row["sample_id"]: {
        "group": row["group"],
        "fastq1": row["fastq1"],
        "fastq2": row["fastq2"]
    }
    for _, row in INCLUDED_SAMPLES.iterrows()
}

# List of included sample IDs
SAMPLES = list(SAMPLE_DICT.keys())

def get_fastq1(wildcards):
    return SAMPLE_DICT[wildcards.sample]["fastq1"]

def get_fastq2(wildcards):
    return SAMPLE_DICT[wildcards.sample]["fastq2"]





include: "rules/fastp.smk"
include: "rules/trinity.smk"
include: "rules/qc.smk"
include: "rules/kallisto.smk"
include: "rules/quantification.smk"
include: "rules/reports.smk"

rule all:
    input:
        expand("results/fastp/{sample}_R1.clean.fastq.gz", sample=SAMPLES),
        expand("results/fastp/{sample}_R2.clean.fastq.gz", sample=SAMPLES),
        expand("results/summary_qc/{sample}_fastp.html", sample=SAMPLES),
        expand("results/summary_qc/{sample}_fastp.json", sample=SAMPLES),
        "results/trinity/trinity_assembly.fasta",
        "results/summary_qc/trinity_metrics_mqc.tsv",
        expand("results/summary_qc/{sample}_bowtie2.log", sample=SAMPLES),
        "results/summary_qc/short_summary.specific.txt",
        "results/rnaquast/short_report.txt",
        expand("results/kallisto/{sample}/abundance.h5", sample=SAMPLES),
        "results/quant/gene_count_matrix.csv",
        "results/quant/transcript_count_matrix.csv",
        "results/summary_qc/pca_plot_mqc.html",
        "results/summary_qc/multiqc_report.html"
