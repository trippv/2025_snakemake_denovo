rule trinity_assembly:
    input:
        r1=expand("results/fastp/{sample}_R1.clean.fastq.gz", sample=SAMPLES),
        r2=expand("results/fastp/{sample}_R2.clean.fastq.gz", sample=SAMPLES)
    output:
        fasta="results/trinity/trinity_assembly.fasta",
        gene_map="results/trinity/trinity_assembly.fasta.gene_trans_map"
    conda:
        "../envs/trinity.yaml"
    threads: config["trinity_threads"]
    resources:
        mem_mb=config["trinity_mem_mb"]
    log:
        "logs/trinity_assembly.log"
    shell:
        """
        Trinity --seqType fq \
                --left {input.r1} \
                --right {input.r2} \
                --CPU {threads} \
                --max_memory {resources.mem_mb}G \
                --bflyCPU 2 \ #para evitar conflictos de paralelismo
                #--full_cleanup \
                --output trinity > {log} 2>&1
        mkdir -p results/trinity
        mv trinity.Trinity.fasta {output.fasta}
        mv trinity.Trinity.fasta.gene_trans_map {output.gene_map}
        """
