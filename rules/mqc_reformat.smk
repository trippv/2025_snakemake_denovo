rule mqc_reformat:
    input:
        abundance = "results/summary_qc/abundance.txt",
        pca = "results/summary_qc/pca.txt",
        volcano = "results/summary_qc/volcano.txt"
    output:
        abundance_yaml = "results/summary_qc/abundance_heatmap_mqc.yaml",
        pca_yaml = "results/summary_qc/pca_mqc.yaml",
        volcano_yaml = "results/summary_qc/volcano_mqc.yaml"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/mqc_reformat_yaml.R"
