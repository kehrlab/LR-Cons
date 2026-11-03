rule create_consensus_sequences:
    input:
        molephase_dir = config["simdata_dir"] + "/{coverage}/reads/{technology}_{accuracy}/molephase",
        fa_list = config["simdata_dir"] + "/{coverage}/{mode}/fasta_list_{technology}_{accuracy}.tsv",
        binary = config["binary_dir"] + "/create_consensus"
    output:
        consensus_dir = directory(config["simdata_dir"] + "/{coverage, [0-9]+}/{mode}/consensus_seqs/{type}/{technology}_{accuracy}/{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}/{method}"),
        consensus_list = config["simdata_dir"] + "/{coverage, [0-9]+}/{mode}/consensus_seqs_{method}_{technology}_{accuracy}_{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}_{type}.tsv",
        stats_out = config["simdata_dir"] + "/{coverage, [0-9]+}/{mode}/consensus_stats_{method}_{technology}_{accuracy}_{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}_{type}.tsv"
    log: "logs/create_consensus_{method}_{coverage}_{technology}_{accuracy}_{p}_{mode}_{type}.log"
    benchmark: "benchmarks/create_consensus_{method}_{coverage}_{technology}_{accuracy}_{p}_{mode}_{type}.txt"
    threads: 8
    resources:
        mem_mb = 8192,
        cpus_per_task = 8
    params:
        phasing_dir = lambda w: config["simdata_dir"] + "/" + w.coverage + "/reads/" + w.technology + "_" + w.accuracy + "/molephase",
        match = matchParam,
        mismatch = mismatchParam,
        gap = openParam,
        ext = extParam,
        k = getK
    shell:
        "{input.binary} {input.fa_list} {output.consensus_dir} {output.consensus_list} {output.stats_out} {params.phasing_dir} {params.match} {params.mismatch} {params.gap} {params.ext} {wildcards.type} {params.k} {threads} 500 {wildcards.method} &> {log}"

rule determine_consensus_times:
    input:
        fa_list = config["simdata_dir"] + "/{coverage}/partial/fasta_list_{technology}_{accuracy}.tsv",
        binary = config["binary_dir"] + "/create_consensus"
    output:
        consensus_dir = temp(directory(config["simdata_dir"] + "/time_data/{coverage, [0-9]+}/consensus_seqs/{technology}_{accuracy}/{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}/{method}")),
        consensus_list = temp(config["simdata_dir"] + "/time_data/{coverage, [0-9]+}/consensus_seqs_{method}_{technology}_{accuracy}_{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}.tsv"),
        stats_out = config["simdata_dir"] + "/time_data/{coverage, [0-9]+}/consensus_stats_{method}_{technology}_{accuracy}_{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}.tsv"
    log: "logs/determine_consensus_times_{method}_{coverage}_{technology}_{accuracy}_{p}.log"
    benchmark: "benchmarks/determine_consensus_times_{method}_{coverage}_{technology}_{accuracy}_{p}.txt"
    threads: 1
    params:
        match = matchParam,
        mismatch = mismatchParam,
        gap = openParam,
        ext = extParam,
        k = getK
    shell:
        "{input.binary} {input.fa_list} {output.consensus_dir} {output.consensus_list} {output.stats_out} /dev/null {params.match} {params.mismatch} {params.gap} {params.ext} local {params.k} 1 500 {wildcards.method} &> {log}"


rule evaluate_consensus_sequences:
    input:
        consensus_list = config["simdata_dir"] + "/{coverage}/{mode}/consensus_seqs_{method}_{technology}_{accuracy}_{p}_{type}.tsv",
        truth_dir = config["simdata_dir"] + "/allele_sequences",
        binary = config["binary_dir"] + "/evaluate_consensus_seqs"
    output:
        config["results_dir"] + "/simulation/{technology}/{accuracy}/{coverage, [0-9]+}/{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}/evaluation_{method}_{mode}_{type}.tsv"
    threads: 8
    resources:
        mem_mb = 4096,
        cpus_per_task = 8
    log: "logs/eval_consensus_{method}_{coverage}_{technology}_{accuracy}_{p}_{mode}_{type}.log"
    shell:
        '{input.binary} {input.consensus_list} {input.truth_dir} {output} {threads} &> {log}'
