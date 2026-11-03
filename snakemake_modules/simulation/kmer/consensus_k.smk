rule create_consensus_sequence_kmer_k:
    input:
        fa_list = config["simdata_dir"] + "/{coverage}/{mode}/fasta_list_{technology}_{accuracy}.tsv",
        binary = config["binary_dir"] + "/create_consensus"
    output:
        consensus_dir = directory(config["simdata_dir"] + "/{coverage, [0-9]+}/{mode}/consensus_seqs/{type}/{technology}_{accuracy}/{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}/KMEANS_{k}"),
        consensus_list = config["simdata_dir"] + "/{coverage, [0-9]+}/{mode}/consensus_seqs_KMEANS_kmer_{k}_{technology}_{accuracy}_{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}_{type}.tsv",
        stats_out = config["simdata_dir"] + "/{coverage, [0-9]+}/{mode}/consensus_stats_KMEANS_kmer_{k}_{technology}_{accuracy}_{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}_{type}.tsv"
    log: "logs/create_consensus_KMEANS_{k}_{coverage}_{technology}_{accuracy}_{p}_{mode}_{type}.log"
    benchmark: "benchmarks/create_consensus_KMEANS_{k}_{coverage}_{technology}_{accuracy}_{p}_{mode}_{type}.txt"
    threads: 16
    params:
        match = matchParam,
        mismatch = mismatchParam,
        gap = openParam,
        ext = extParam
    shell:
        "{input.binary} {input.fa_list} {output.consensus_dir} {output.consensus_list} {output.stats_out} /dev/null {params.match} {params.mismatch} {params.gap} {params.ext} {wildcards.type} {wildcards.k} {threads} 500 KMEANS &> {log}"

rule evaluate_consensus_sequences_kmer_k:
    input:
        consensus_list = config["simdata_dir"] + "/{coverage}/{mode}/consensus_seqs_KMEANS_kmer_{k}_{technology}_{accuracy}_{p}_{type}.tsv",
        truth_dir = config["simdata_dir"] + "/allele_sequences",
        binary = config["binary_dir"] + "/evaluate_consensus_seqs"
    output:
        config["results_dir"] + "/kmeans_simulation/{technology}/{accuracy}/{coverage, [0-9]+}/{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}/evaluation_KMEANS_{k}_{mode}_{type}.tsv"
    threads: 15
    log: "logs/eval_consensus_KMEANS_{k}_{coverage}_{technology}_{accuracy}_{p}_{mode}_{type}.log"
    shell:
        '{input.binary} {input.consensus_list} {input.truth_dir} {output} {threads} &> {log}'
