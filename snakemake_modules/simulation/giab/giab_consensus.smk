rule create_consensus_sequences_giab_sim:
    input:
        molephase_dir = config["simdata_dir"] + "/giab/{technology}/{coverage}/reads/{technology}_0.99/molephase",
        fa_list = config["simdata_dir"] + "/giab/{technology}/{coverage}/fasta_list_{technology}_{mode}.tsv",
        binary = config["binary_dir"] + "/create_consensus"
    output:
        consensus_dir = directory(config["simdata_dir"] + "/giab/{technology}/{coverage}/consensus_seqs/{mode}/{type}/{method}"),
        consensus_list = config["simdata_dir"] + "/giab/{technology}/{coverage}/consensus_seqs_{method}_{mode}_{type}.tsv",
        stats_out = config["simdata_dir"] + "/giab/{technology}/{coverage}/consensus_stats_{method}_{mode}_{type}.tsv"
    log: "logs/giab_sim_create_consensus_{method}_{technology}_{mode}_{type}_{coverage}.log"
    benchmark: "benchmarks/giab_sim_create_consensus_{method}_{technology}_{mode}_{type}_{coverage}.txt"
    threads: 8
    resources:
        mem_mb = 12000,
        cpus_per_task = 8
    priority: 4
    params:
        #iphasing_dir = lambda w: config["simdata_dir"] + "/giab/" + w.technology + "/" + w.coverage + "/reads/" + w.technology + "_0.99/molephase",
        k = getK,
        match = 1,
        mismatch = -1,
        gap = -1,
        ext = -1
    shell:
        "{input.binary} {input.fa_list} {output.consensus_dir} {output.consensus_list} {output.stats_out} {input.molephase_dir} {params.match} {params.mismatch} {params.gap} {params.ext} {wildcards.type} {params.k} {threads} 500 {wildcards.method} &> {log}"

rule evaluate_consensus_sequences_giab_sim:
    input:
        consensus_list = config["simdata_dir"] + "/giab/{technology}/{coverage}/consensus_seqs_{method}_{mode}_{type}.tsv",
        truth_dir = config["simdata_dir"] + "/giab/allele_sequences",
        binary = config["binary_dir"] + "/evaluate_consensus_seqs"
    output:
        config["results_dir"] + "/giab_simulation/{technology}/{coverage}/evaluation_{method}_{mode}_{type}.tsv"
    threads: 8
    priority: 3
    resources:
        mem_mb = 8000,
        cpus_per_task = 8
    log: "logs/giab_sim_eval_consensus_{method}_{technology}_{mode}_{type}_{coverage}.log"
    shell:
        '{input.binary} {input.consensus_list} {input.truth_dir} {output} {threads} &> {log}'
