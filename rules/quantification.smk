rule build_kallisto_sample_table:
    input:
        abundances=expand("results/kallisto/{sample}/abundance.tsv", sample=SAMPLES)
    output:
        table="results/quant/kallisto_samples_table.txt"
    run:
        with open(output.table, 'w') as f:
            for sample in SAMPLES:
                f.write(f"{sample}\tresults/kallisto/{sample}/abundance.tsv\n")

rule generate_metadata:
    input:
        sample_file=config["samples_file"]
    output:
        metadata="data/metadata.tsv"
    run:
        import pandas as pd
        samples = pd.read_csv(input.sample_file, sep="\t")
        samples = samples[samples["include"] == 1][["sample_id", "group"]]
        samples.to_csv(output.metadata, sep="\t", index=False,
                       header=["sample", "group"])

rule kallisto_to_matrix:
    input:
        table="results/quant/kallisto_samples_table.txt",
        gene_trans_map="results/trinity/trinity_assembly.fasta.gene_trans_map"
    output:
        gene_matrix="results/quant/gene_count_matrix.csv",
        transcript_matrix="results/quant/transcript_count_matrix.csv"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/kallisto_to_matrix.R"
