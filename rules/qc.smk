rule busco:
    input:
        fasta="results/trinity/trinity_assembly.fasta"
    output:
        summary="results/busco/short_summary.specific.txt"
    conda:
        "../envs/busco.yaml"
    threads: config["busco_threads"]
    resources:
        mem_mb=config["busco_mem_mb"]
    params:
        lineage=config["busco_lineage"],
        out_dir="results/busco"
    log:
        "logs/busco/busco.log"
    shell:
        """
        busco -i {input.fasta} \
              -o busco_output \
              -m transcriptome \
              -l {params.lineage} \
              -c {threads} \
              -f \
              --out_path {params.out_dir} > {log} 2>&1
        mv {params.out_dir}/busco_output/short_summary*.txt {output.summary}
        """

rule rnaquast:
    input:
        fasta="results/trinity/trinity_assembly.fasta"
    output:
        report="results/rnaquast/short_report.txt"
    conda:
        "../envs/rnaquast.yaml"
    threads: config.get("rnaquast_threads", 2)
    resources:
        mem_mb=config.get("rnaquast_mem_mb", 2000)
    params:
        out_dir="results/rnaquast/"
    log:
        "logs/rnaquast/rnaquast.log"
    shell:
        """
        rnaQUAST.py --transcripts {input.fasta} \
                    --output_dir {params.out_dir} \
                    --threads {threads} \
                    --no_plots > {log} 2>&1
        """
rule trinity_stats:
    input:
        fasta = "results/trinity/trinity_assembly.fasta"
    output:
        stats = "results/trinity/trinity_stats.txt"
    conda:
        "../envs/trinity.yaml"
    log:
        "logs/trinity_stats.log"
    shell:
        """
        TrinityStats.pl {input.fasta} > {output.stats} 2> {log}
        """

rule parse_trinity_stats:
    input:
        stats = "results/trinity/trinity_stats.txt"
    output:
        tsv = "results/trinity/trinity_metrics_mqc.tsv"
    run:
        import re

        with open(input.stats, 'r') as f:
            content = f.read()

        def quick_find(pattern, text):
            match = re.search(pattern, text)
            return match.group(1).strip() if match else "0"

        # 1. Extraer datos generales (Sección superior)
        total_genes = quick_find(r"Total trinity 'genes':\s+(\d+)", content)
        total_trans = quick_find(r"Total trinity transcripts:\s+(\d+)", content)
        gc_content  = quick_find(r"Percent GC:\s+([\d\.]+)", content)

        # 2. Extraer métricas de ALL transcript contigs (Primera sección de stats)
        # Usamos un re.search limitado a la primera parte del archivo para evitar confundir con 'Longest Isoform'
        all_stats_part = content.split("Stats based on ONLY LONGEST ISOFORM")[0]
        
        n10         = quick_find(r"Contig N10:\s+(\d+)", all_stats_part)
        n50         = quick_find(r"Contig N50:\s+(\d+)", all_stats_part)
        median_len  = quick_find(r"Median contig length:\s+(\d+)", all_stats_part)
        avg_len     = quick_find(r"Average contig:\s+([\d\.]+)", all_stats_part)
        total_bases = quick_find(r"Total assembled bases:\s+(\d+)", all_stats_part)

        # 3. Crear el archivo TSV para MultiQC
        with open(output.tsv, 'w') as out:
            out.write("# id: 'trinity_stats'\n")
            out.write("# section_name: 'Estadísticas del Ensamblaje (Trinity)'\n")
            out.write("# plot_type: 'table'\n")
            out.write("# pconfig:\n")
            out.write("#    namespace: 'Trinity Metrics'\n")
            
            # Cabecera con más detalles
            header = ["Sample", "Genes", "Transcripts", "N10", "N50", "Median_Len", "Avg_Len", "Total_Bases", "%GC"]
            out.write("\t".join(header) + "\n")
            
            # Valores
            values = ["Trinity_Assembly", total_genes, total_trans, n10, n50, median_len, avg_len, total_bases, gc_content]
            out.write("\t".join(values) + "\n")


# Bowtie para mapear las lacturas al ensamblaje de Trinity
# 1. Crear el índice de Bowtie2 a partir del ensamblaje de Trinity
rule bowtie2_build_index:
    input:
        fasta = "results/trinity/trinity_assembly.fasta"
    output:
        # Bowtie2 genera varios archivos .bt2. Usamos un "touch" o el primer archivo como referencia.
        directory("results/mapping/bowtie2_index"),
        index_base = "results/mapping/bowtie2_index/trinity_idx.1.bt2"
    log:
        "logs/bowtie2_index.log"
    conda:
        "../envs/bowtie.yaml"
    params:
        prefix = "results/mapping/bowtie2_index/trinity_idx"
    shell:
        """
        mkdir -p results/mapping/bowtie2_index
        bowtie2-build {input.fasta} {params.prefix} > {log} 2>&1
        """

# 2. Alinear las lecturas procesadas (fastp) contra el ensamblaje
rule bowtie2_align_to_assembly:
    input:
        r1 = "results/fastp/{sample}_R1.clean.fastq.gz",
        r2 = "results/fastp/{sample}_R2.clean.fastq.gz",
        index = "results/mapping/bowtie2_index/trinity_idx.1.bt2"
    output:
        bam = "results/mapping/{sample}_to_assembly.bam",
        log = "results/summary_qc/{sample}_bowtie2.log"
    params:
        prefix = "results/mapping/bowtie2_index/trinity_idx"
    conda:
        "../envs/bowtie.yaml"
    threads: 8
    shell:
        """
        bowtie2 -p {threads} -x {params.prefix} \
            -1 {input.r1} -2 {input.r2} \
            --no-unal --sensitive --no-mixed --no-discordant \
            2> {output.log} | samtools view -bS - > {output.bam}
        """