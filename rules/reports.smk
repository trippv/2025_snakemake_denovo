rule pca_html:
    input:
        counts="results/quant/gene_count_matrix.csv",
        metadata="data/metadata.tsv",
        rmd="scripts/pca_plot.Rmd"
    output:
        html="results/summary_qc/pca_plot_mqc.html"
    conda:
        "../envs/r.yaml"
    shell:
        """
        MAIN=$(pwd)
        Rscript -e "rmarkdown::render('{input.rmd}', params=list(main='$MAIN', counts='{input.counts}', metadata='{input.metadata}'))"
        mv scripts/pca_plot.html {output.html}
        """

rule multiqc:
    input:
        fastp=expand("results/summary_qc/{sample}_fastp.json", sample=SAMPLES),
        busco="results/busco/short_summary.specific.txt"
    output:
        report="results/summary_qc/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    params:
        comment=config["multiqc_comment"]
    shell:
        """
        multiqc results/summary_qc -f \
            --config config/multiqc_config.yaml \
            --template mi_multiqc \
            --comment "{params.comment}"
        mv multiqc_report.html {output.report}
        mv multiqc_data results/summary_qc/
        """
