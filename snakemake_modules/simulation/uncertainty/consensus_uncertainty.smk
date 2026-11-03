rule create_consensus_sequences_uncertainty:
    input:
        molephase_dir = config["simdata_dir"] + "/15/reads/{technology}_0.99/molephase",
        fa_list = config["simdata_dir"] + "/uncertainty_{margin}_{jitter}/{mode}/fasta_list_{technology}.tsv",
        binary = config["binary_dir"] + "/create_consensus"
    output:
        consensus_dir = directory(config["simdata_dir"] + "/uncertainty_{margin}_{jitter}/{mode}/consensus_seqs/{type}/{technology}/{method}"),
        consensus_list = config["simdata_dir"] + "/uncertainty_{margin}_{jitter}/{mode}/consensus_seqs_{method}_{technology}_{type}.tsv",
        stats_out = config["simdata_dir"] + "/uncertainty_{margin}_{jitter}/{mode}/consensus_stats_{method}_{technology}_{type}.tsv"
    log: "logs/create_consensus_{method}_{technology}_{mode}_{type}_{margin}_{jitter}.log"
    benchmark: "benchmarks/create_consensus_{method}_{technology}_{mode}_{type}_{margin}_{jitter}.txt"
    threads: 8
    resources:
        mem_mb = 12000
    priority: 2
    params:
        phasing_dir = lambda w: config["simdata_dir"] + "/15/reads/" + w.technology + "_0.99/molephase",
        k = getK,
        match = 1,
        mismatch = -1,
        gap = -1,
        ext = -1
    shell:
        "{input.binary} {input.fa_list} {output.consensus_dir} {output.consensus_list} {output.stats_out} {params.phasing_dir} {params.match} {params.mismatch} {params.gap} {params.ext} {wildcards.type} {params.k} {threads} {wildcards.margin} {wildcards.method} &> {log}"

rule evaluate_consensus_sequences_uncertainty:
    input:
        consensus_list = config["simdata_dir"] + "/uncertainty_{margin}_{jitter}/{mode}/consensus_seqs_{method}_{technology}_{type}.tsv",
        truth_dir = config["simdata_dir"] + "/allele_sequences",
        binary = config["binary_dir"] + "/evaluate_consensus_seqs"
    output:
        config["results_dir"] + "/uncertainty_simulation/{technology}/{margin}/{jitter}/evaluation_{method}_{mode}_{type}.tsv"
    threads: 8
    priority: 1
    resources:
        mem_mb = 8000
    log: "logs/eval_consensus_{method}_{technology}_{mode}_{type}_{margin}_{jitter}.log"
    shell:
        '{input.binary} {input.consensus_list} {input.truth_dir} {output} {threads} &> {log}'
