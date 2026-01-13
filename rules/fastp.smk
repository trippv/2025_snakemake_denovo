rule fastp:
    input:
        r1=get_fastq1,
        r2=get_fastq2
    output:
        r1_clean="results/fastp/{sample}_R1.clean.fastq.gz",
        r2_clean="results/fastp/{sample}_R2.clean.fastq.gz",
        json="results/summary_qc/{sample}_fastp.json",
        html="results/summary_qc/{sample}_fastp.html"
    conda:
        "../envs/fastp.yaml"
    threads: config.get("fastp_threads", 2)
    resources:
        mem_mb=config.get("fastp_mem_mb", 2000)
    shell:
        """
        fastp --in1 {input.r1} --in2 {input.r2} \
              --out1 {output.r1_clean} --out2 {output.r2_clean} \
              --json {output.json} --html {output.html} \
              --thread {threads} \
              --detect_adapter_for_pe
        """
