rule expr_diff:
    input:
        counts = "results/quant/gene_count_matrix.csv",
        metadata = "data/metadata.tsv"
    output:
        de_dir = directory("results/Exp_Diff"),
        dist_matrix="results/summary_qc/distance_matrix.txt",
        pca="results/summary_qc/pca.txt",
        abundance="results/summary_qc/abundance.txt",
        volcano="results/summary_qc/volcano.txt"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/ExprDiff.R"
