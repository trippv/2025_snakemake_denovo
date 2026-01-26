rule busco:
    input:
        fasta="results/trinity/trinity_assembly.fasta"
    output:
        summary="results/summary_qc/short_summary.specific.txt"
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
    shell:
        """
        TrinityStats.pl {input.fasta} > {output.stats}
        """

#rule parse_trinity_stats:
#    input:
#        stats = "results/trinity/trinity_stats.txt"
#    output:
#        tsv = "results/summary_qc/trinity_metrics_mqc.tsv"
#    run:
#        import re
#
#        with open(input.stats, 'r') as f:
#            content = f.read()
#
#        def quick_find(pattern, text):
#            match = re.search(pattern, text)
#            return match.group(1).strip() if match else "0"
#
#        # 1. Extraer datos generales (Sección superior)
#        total_genes = quick_find(r"Total trinity 'genes':\s+(\d+)", content)
#        total_trans = quick_find(r"Total trinity transcripts:\s+(\d+)", content)
#        gc_content  = quick_find(r"Percent GC:\s+([\d\.]+)", content)
#
#        # 2. Extraer métricas de ALL transcript contigs (Primera sección de stats)
#        # Usamos un re.search limitado a la primera parte del archivo para evitar confundir con 'Longest Isoform'
#        all_stats_part = content.split("Stats based on ONLY LONGEST ISOFORM")[0]
#        
#        n10         = quick_find(r"Contig N10:\s+(\d+)", all_stats_part)
#        n50         = quick_find(r"Contig N50:\s+(\d+)", all_stats_part)
#        median_len  = quick_find(r"Median contig length:\s+(\d+)", all_stats_part)
#        avg_len     = quick_find(r"Average contig:\s+([\d\.]+)", all_stats_part)
#        total_bases = quick_find(r"Total assembled bases:\s+(\d+)", all_stats_part)
#
#        # 3. Crear el archivo TSV para MultiQC
#        with open(output.tsv, 'w') as out:
#            out.write("# id: 'trinity_stats'\n")
#            out.write("# section_name: 'Estadísticas del Ensamblaje (Trinity)'\n")
#            out.write("# plot_type: 'table'\n")
#            out.write("# pconfig:\n")
#            out.write("#    namespace: 'Trinity Metrics'\n")
#            
#            # Cabecera con más detalles
#            header = ["Sample", "Genes", "Transcripts", "N10", "N50", "Median_Len", "Avg_Len", "Total_Bases", "%GC"]
#            out.write("\t".join(header) + "\n")
#            
#            # Valores
#            values = ["Trinity_Assembly", total_genes, total_trans, n10, n50, median_len, avg_len, total_bases, gc_content]
#            out.write("\t".join(values) + "\n")


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

rule parse_assembly_quality:
    input:
        stats = "results/trinity/trinity_stats.txt",
        fasta = "results/trinity/trinity_assembly.fasta",
        gff3  = "results/transdecoder/trinity_assembly.fasta.transdecoder.gff3"
    output:
        tsv = "results/summary_qc/trinity_metrics_mqc.tsv"
    run:
        import re

        # --- 1. PROCESAR TRINITY STATS ---
        with open(input.stats, 'r') as f:
            content = f.read()

        def quick_find(pattern, text):
            match = re.search(pattern, text)
            return match.group(1).strip() if match else "0"

        total_genes = quick_find(r"Total trinity 'genes':\s+(\d+)", content)
        total_trans = quick_find(r"Total trinity transcripts:\s+(\d+)", content)
        gc_content  = quick_find(r"Percent GC:\s+([\d\.]+)", content)

        all_stats_part = content.split("Stats based on ONLY LONGEST ISOFORM")[0]
        n50         = quick_find(r"Contig N50:\s+(\d+)", all_stats_part)
        avg_len     = quick_find(r"Average contig:\s+([\d\.]+)", all_stats_part)
        total_bases = quick_find(r"Total assembled bases:\s+(\d+)", all_stats_part)

        # --- 2. CALCULAR POTENCIAL CODIFICANTE (TRANSDECODER) ---
        ids_con_orf = set()
        with open(input.gff3, 'r') as f:
            for line in f:
                if not line.startswith("#") and "\tmRNA\t" in line:
                    # El ID del transcripto original en TransDecoder está en la columna 1
                    ids_con_orf.add(line.split("\t")[0])
        
        n_orfs = len(ids_con_orf)
        # Usamos total_trans convertido a float para el cálculo
        percent_orf = (n_orfs / float(total_trans) * 100) if float(total_trans) > 0 else 0

        # --- 3. ESCRIBIR TSV UNIFICADO ---
        with open(output.tsv, 'w') as out:
            out.write("# id: 'assembly_qc_stats'\n")
            out.write("# section_name: 'Calidad Integral del Ensamblaje'\n")
            out.write("# plot_type: 'table'\n")
            out.write("# pconfig:\n")
            out.write("#     namespace: 'Assembly Metrics'\n")
            
            # Cabecera expandida con métricas de TransDecoder
            header = ["Sample", "Genes", "Transcripts", "N50", "Avg_Len", "Total_MB", "%GC", "Trans_with_ORF", "%_with_ORF"]
            out.write("\t".join(header) + "\n")
            
            # Convertir bases a Megabases para que sea más legible (opcional)
            mb_assembled = f"{int(total_bases)/1e6:.2f}"
            
            # Valores combinados
            values = [
                "Trinity_Assembly", 
                total_genes, 
                total_trans, 
                n50, 
                avg_len, 
                mb_assembled, 
                gc_content, 
                str(n_orfs), 
                f"{percent_orf:.2f}"
            ]
            out.write("\t".join(values) + "\n")