rule create_consensus_sequences_real:
    input:
        directory(config["realdata_dir"] + "/reads/{technology}_{coverage}/molephase"),
        fa_list = config["realdata_dir"] + "/{mode}/fasta_list_{technology}_{coverage}_{margin}.tsv",
        binary = config["binary_dir"] + "/create_consensus"
    output:
        consensus_dir = directory(config["realdata_dir"] + "/{mode}/{type}/consensus_seqs/{technology}/{method}/{coverage}_{margin}"),
        consensus_list = config["realdata_dir"] + "/{mode}/{type}/consensus_seqs_{method}_{technology}_{coverage}_{margin}.tsv",
        stats_out = config["realdata_dir"] + "/{mode}/{type}/consensus_stats_{method}_{technology}_{coverage}_{margin}.tsv"
    log: "logs/create_consensus_{mode}_{type}_{method}_{technology}_{coverage}_{margin}.log"
    benchmark: "benchmarks/create_consensus_{mode}_{type}_{method}_{technology}_{coverage}_{margin}.txt"
    threads: 8
    priority: 2
    resources:
        mem_mb = 20000
    params:
        phasing_dir = lambda w: config["realdata_dir"] + "/reads/" + w.technology + "_" + w.coverage + "/molephase",
        match = "1",
        mismatch = "-1",
        gap = "-1",
        ext = "-1",
        m = 500,
        k = getK
    shell:
        "{input.binary} {input.fa_list} {output.consensus_dir} {output.consensus_list} {output.stats_out} {params.phasing_dir} {params.match} {params.mismatch} {params.gap} {params.ext} {wildcards.type} {params.k} {threads} {params.m} {wildcards.method} &> {log}"

rule time_benchmark_consensus_real:
    input:
        fa_list = config["realdata_dir"] + "/{mode}/fasta_list_{technology}_{coverage}_{margin}.tsv",
        binary = config["binary_dir"] + "/create_consensus"
    output:
        consensus_dir = directory("data/time_benchmark/{mode}/{type}/consensus_seqs/{technology}/{method}/{coverage}_{margin}"),
        consensus_list = "data/time_benchmark/{mode}/{type}/consensus_seqs_{method}_{technology}_{coverage}_{margin}.tsv",
        stats_out = "data/time_benchmark/{mode}/{type}/consensus_stats_{method}_{technology}_{coverage}_{margin}.tsv"
    log: "logs/time_create_consensus_{mode}_{type}_{method}_{technology}_{coverage}_{margin}.log"
    benchmark: "benchmarks/time_create_consensus_{mode}_{type}_{method}_{technology}_{coverage}_{margin}.txt"
    threads: 1
    priority: 2
    resources:
        mem_mb = 20000
    params:
        match = "1",
        mismatch = "-1",
        gap = "-1",
        ext = "-1",
        m = 500,
        k = getK
    shell:
        "{input.binary} {input.fa_list} {output.consensus_dir} {output.consensus_list} {output.stats_out} /dev/null {params.match} {params.mismatch} {params.gap} {params.ext} {wildcards.type} {params.k} {threads} {params.m} {wildcards.method} &> {log}"



rule evaluate_consensus_sequences_real:
    input:
        consensus_list = config["realdata_dir"] + "/{extraction}/{mode}/consensus_seqs_{method}_{technology}_{coverage}_{margin}.tsv",
        truth_dir = config["asm_dir"] + '/variant_sequences_{margin}',
        binary = config["binary_dir"] + "/evaluate_consensus_seqs"
    output:
        config["results_dir"] + "/real/{technology}/{coverage}/evaluation_{method}_{extraction}_{mode}_{margin, [0-9]+}.tsv"
    threads: 8
    priority: 1
    resources:
        mem_mb = 8192
    log: "logs/eval_consensus_{method}_{coverage}_{technology}_{extraction}_{mode}_{margin}.log"
    shell:
        '{input.binary} {input.consensus_list} {input.truth_dir} {output} {threads} &> {log}'
