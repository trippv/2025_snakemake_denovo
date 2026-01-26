rule transdecoder_longorfs:
    input:
        fasta = "results/trinity/trinity_assembly.fasta"
    output:
        # El directorio que ya viste que se crea
        td_dir = directory("results/transdecoder/trinity_assembly.fasta.transdecoder_dir"),
        # El archivo pep que está dentro (podemos referenciarlo directamente)
        pep_long = "results/transdecoder/trinity_assembly.fasta.transdecoder_dir/longest_orfs.pep"
    log:
        "logs/transdecoder/longorfs.log"
#    conda:
#        "../envs/transdecoder.yaml"
    shell:
        """
        TransDecoder.LongOrfs -t {input.fasta} --output_dir results/transdecoder > {log} 2>&1
        """

rule transdecoder_predict:
    input:
        fasta = "results/trinity/trinity_assembly.fasta",
        td_dir = "results/transdecoder/trinity_assembly.fasta.transdecoder_dir"
    output:
        # Los archivos finales que aparecen en tu listado
        pep = "results/transdecoder/trinity_assembly.fasta.transdecoder.pep",
        gff3 = "results/transdecoder/trinity_assembly.fasta.transdecoder.gff3"
    log:
        "logs/transdecoder/predict.log"
#    conda:
#        "../envs/transdecoder.yaml"
    shell:
        """
        TransDecoder.Predict -t {input.fasta} --output_dir results/transdecoder > {log} 2>&1
        """
