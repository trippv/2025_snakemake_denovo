rule transdecoder_longorfs:
    input:
        fasta = "results/trinity/trinity_assembly.fasta"
    output:
        # El archivo clave que Predict necesita
        pep_long = "results/transdecoder/longest_orfs.pep"
    log:
        "logs/transdecoder/longorfs.log"
    conda:
        "../envs/transdecoder.yaml"
    shell:
        """
        # 1. Crear directorio si no existe
        mkdir -p results/transdecoder/

        # 2. Ejecutar LongOrfs
        TransDecoder.LongOrfs -t {input.fasta} --output_dir results/transdecoder > {log} 2>&1
        """

rule transdecoder_predict:
    input:
        fasta = "results/trinity/trinity_assembly.fasta",
        pep_long = "results/transdecoder/longest_orfs.pep"
    output:
        pep = "results/transdecoder/trinity_assembly.fasta.transdecoder.pep",
        gff3 = "results/transdecoder/trinity_assembly.fasta.transdecoder.gff3",
        cds = "results/transdecoder/trinity_assembly.fasta.transdecoder.cds",
        bed = "results/transdecoder/trinity_assembly.fasta.transdecoder.bed"
    log:
        "logs/transdecoder/predict.log"
    conda:
        "../envs/transdecoder.yaml"
    shell:
        """
        # 1. Ejecutar Predict usando el mismo output_dir que LongOrfs
        TransDecoder.Predict -t {input.fasta} --output_dir results/transdecoder >> {log} 2>&1

        # 2. TransDecoder dejó los archivos en la raíz (donde se ejecuta snakemake).
        # Los movemos manualmente a la carpeta results/transdecoder/
        
        NAME=$(basename {input.fasta})

        mv ${{NAME}}.transdecoder.pep {output.pep} >> {log} 2>&1
        mv ${{NAME}}.transdecoder.gff3 {output.gff3} >> {log} 2>&1
        mv ${{NAME}}.transdecoder.cds {output.cds} >> {log} 2>&1
        mv ${{NAME}}.transdecoder.bed {output.bed} >> {log} 2>&1
        """