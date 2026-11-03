rule create_consensus_sequences_real_old:
    input:
        directory(config["realdata_dir"] + "/old/reads/{technology}_{coverage}/molephase"),
        fa_list = config["realdata_dir"] + "/old/{mode}/fasta_list_{technology}_{coverage}_{margin}.tsv",
        binary = config["binary_dir"] + "/create_consensus"
    output:
        consensus_dir = directory(config["realdata_dir"] + "/old/{mode}/{type}/consensus_seqs/{technology}/{method}/{coverage}_{margin}"),
        consensus_list = config["realdata_dir"] + "/old/{mode}/{type}/consensus_seqs_{method}_{technology}_{coverage}_{margin}.tsv",
        stats_out = config["realdata_dir"] + "/old/{mode}/{type}/consensus_stats_{method}_{technology}_{coverage}_{margin}.tsv"
    log: "logs/create_consensus_{mode}_{type}_{method}_{technology}_{coverage}_{margin}_old.log"
    benchmark: "benchmarks/create_consensus_{mode}_{type}_{method}_{technology}_{coverage}_{margin}_old.txt"
    threads: 8
    priority: 2
    resources:
        mem_mb = 20000
    params:
        phasing_dir = lambda w: config["realdata_dir"] + "/old/reads/" + w.technology + "_" + w.coverage + "/molephase",
        match = "1",
        mismatch = "-1",
        gap = "-1",
        ext = "-1",
        m = 500,
        k = getK
    shell:
        "{input.binary} {input.fa_list} {output.consensus_dir} {output.consensus_list} {output.stats_out} {params.phasing_dir} {params.match} {params.mismatch} {params.gap} {params.ext} {wildcards.type} {params.k} {threads} {params.m} {wildcards.method} &> {log}"

rule evaluate_consensus_sequences_real_old:
    input:
        consensus_list = config["realdata_dir"] + "/old/{extraction}/{mode}/consensus_seqs_{method}_{technology}_{coverage}_{margin}.tsv",
        truth_dir = config["asm_dir"] + '/variant_sequences_{margin}',
        binary = config["binary_dir"] + "/evaluate_consensus_seqs"
    output:
        config["results_dir"] + "/real/{technology}_old/{coverage}/evaluation_{method}_{extraction}_{mode}_{margin, [0-9]+}.tsv"
    threads: 8
    priority: 1
    resources:
        mem_mb = 8192
    log: "logs/eval_consensus_{method}_{coverage}_{technology}_{extraction}_{mode}_{margin}_old.log"
    shell:
        '{input.binary} {input.consensus_list} {input.truth_dir} {output} {threads} &> {log}'
