rule kallisto_index:
    input:
        fasta="results/trinity/trinity_assembly.fasta"
    output:
        index="results/kallisto/trinity_index.idx"
    conda:
        "../envs/kallisto.yaml"
    resources:
        mem_mb=config["kallisto_mem_mb"]
    shell:
        "kallisto index -i {output.index} {input.fasta}"

rule kallisto_quant:
    input:
        r1="results/fastp/{sample}_R1.clean.fastq.gz",
        r2="results/fastp/{sample}_R2.clean.fastq.gz",
        index="results/kallisto/trinity_index.idx"
    output:
        abundance="results/kallisto/{sample}/abundance.h5",
        tsv="results/kallisto/{sample}/abundance.tsv"
    conda:
        "../envs/kallisto.yaml"
    threads: config["kallisto_threads"]
    resources:
        mem_mb=config["kallisto_mem_mb"]
    log:
        "results/summary_qc/{sample}.log"
    shell:
        """
        kallisto quant -i {input.index} \
                       -o results/kallisto/{wildcards.sample} \
                       -t {threads} \
                       {input.r1} {input.r2} > {log} 2>&1
        """
