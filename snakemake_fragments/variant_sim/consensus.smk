ruleorder: evaluate_consensus_sequences_global > evaluate_consensus_sequences_semiglobal > evaluate_consensus_sequences_spanning > evaluate_consensus_sequences

rule create_consensus_sequences:
    input:
        expand(config["data_dir"] + "/{{coverage}}/reads/{{technology}}_{{accuracy}}/molephase/{region}.hap1.fq", region = define_phasing_regions()),
        expand(config["data_dir"] + "/{{coverage}}/reads/{{technology}}_{{accuracy}}/molephase/{region}.hap2.fq", region = define_phasing_regions()),
        expand(config["data_dir"] + "/{{coverage}}/reads/{{technology}}_{{accuracy}}/molephase/{region}.unphased.fq", region = define_phasing_regions()),
        fa_list = config["data_dir"] + "/{coverage}/fasta_list_{technology}_{accuracy}.tsv",
        binary = config["binary_dir"] + "create_consensus"
    output:
        consensus_dir = directory(config["data_dir"] + "/{coverage, [0-9]+}/consensus_seqs/{technology}_{accuracy}/{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}/{method}"),
        consensus_list = config["data_dir"] + "/{coverage, [0-9]+}/consensus_seqs_{method}_{technology}_{accuracy}_{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}.tsv",
        stats_out = config["data_dir"] + "/{coverage, [0-9]+}/consensus_stats_{method}_{technology}_{accuracy}_{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}.tsv"
    log: "logs/create_consensus_{method}_{coverage}_{technology}_{accuracy}_{p}.log"
    benchmark: "benchmarks/create_consensus_{method}_{coverage}_{technology}_{accuracy}_{p}.txt"
    threads: 16
    params:
        phasing_dir = lambda w: config["data_dir"] + "/" + w.coverage + "/reads/" + w.technology + "_" + w.accuracy + "/molephase",
        match = matchParam,
        mismatch = mismatchParam,
        gap = openParam,
        ext = extParam
    shell:
        "{input.binary} {input.fa_list} {output.consensus_dir} {output.consensus_list} {output.stats_out} {params.phasing_dir} {params.match} {params.mismatch} {params.gap} {params.ext} local {wildcards.method} &> {log}"


rule create_consensus_sequences_semiglobal:
    input:
        expand(config["data_dir"] + "/{{coverage}}/reads/{{technology}}_{{accuracy}}/molephase/{region}.hap1.fq", region = define_phasing_regions()),
        expand(config["data_dir"] + "/{{coverage}}/reads/{{technology}}_{{accuracy}}/molephase/{region}.hap2.fq", region = define_phasing_regions()),
        expand(config["data_dir"] + "/{{coverage}}/reads/{{technology}}_{{accuracy}}/molephase/{region}.unphased.fq", region = define_phasing_regions()),
        fa_list = config["data_dir"] + "/{coverage}/fasta_list_{technology}_{accuracy}.tsv",
        binary = config["binary_dir"] + "create_consensus"
    output:
        consensus_dir = directory(config["data_dir"] + "/{coverage, [0-9]+}/semiglobal/consensus_seqs/{technology}_{accuracy}/{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}/{method}"),
        consensus_list = config["data_dir"] + "/{coverage, [0-9]+}/semiglobal/consensus_seqs_{method}_{technology}_{accuracy}_{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}.tsv",
        stats_out = config["data_dir"] + "/{coverage, [0-9]+}/semiglobal/consensus_stats_{method}_{technology}_{accuracy}_{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}.tsv"
    log: "logs/create_consensus_semiglobal_{method}_{coverage}_{technology}_{accuracy}_{p}.log"
    benchmark: "benchmarks/create_consensus_semiglobal_{method}_{coverage}_{technology}_{accuracy}_{p}.txt"
    threads: 16
    params:
        phasing_dir = lambda w: config["data_dir"] + "/" + w.coverage + "/reads/" + w.technology + "_" + w.accuracy + "/molephase",
        match = matchParam,
        mismatch = mismatchParam,
        gap = openParam,
        ext = extParam
    shell:
        "{input.binary} {input.fa_list} {output.consensus_dir} {output.consensus_list} {output.stats_out} {params.phasing_dir} {params.match} {params.mismatch} {params.gap} {params.ext} semiglobal {wildcards.method} &> {log}"


rule create_consensus_sequences_spanning:
    input:
        expand(config["data_dir"] + "/{{coverage}}/reads/{{technology}}_{{accuracy}}/molephase/{region}.hap1.fq", region = define_phasing_regions()),
        expand(config["data_dir"] + "/{{coverage}}/reads/{{technology}}_{{accuracy}}/molephase/{region}.hap2.fq", region = define_phasing_regions()),
        expand(config["data_dir"] + "/{{coverage}}/reads/{{technology}}_{{accuracy}}/molephase/{region}.unphased.fq", region = define_phasing_regions()),
        fa_list = config["data_dir"] + "/{coverage}/global/fasta_list_{technology}_{accuracy}.tsv",
        binary = config["binary_dir"] + "create_consensus"
    output:
        consensus_dir = directory(config["data_dir"] + "/{coverage, [0-9]+}/spanning/consensus_seqs/{technology}_{accuracy}/{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}/{method}"),
        consensus_list = config["data_dir"] + "/{coverage, [0-9]+}/spanning/consensus_seqs_{method}_{technology}_{accuracy}_{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}.tsv",
        stats_out = config["data_dir"] + "/{coverage, [0-9]+}/spanning/consensus_stats_{method}_{technology}_{accuracy}_{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}.tsv"
    log: "logs/create_consensus_spanning_{method}_{coverage}_{technology}_{accuracy}_{p}.log"
    benchmark: "benchmarks/create_consensus_spanning_{method}_{coverage}_{technology}_{accuracy}_{p}.txt"
    threads: 16
    params:
        phasing_dir = lambda w: config["data_dir"] + "/" + w.coverage + "/reads/" + w.technology + "_" + w.accuracy + "/molephase",
        match = matchParam,
        mismatch = mismatchParam,
        gap = openParam,
        ext = extParam
    shell:
        "{input.binary} {input.fa_list} {output.consensus_dir} {output.consensus_list} {output.stats_out} {params.phasing_dir} {params.match} {params.mismatch} {params.gap} {params.ext} local {wildcards.method} &> {log}"

rule create_consensus_sequences_global:
    input:
        expand(config["data_dir"] + "/{{coverage}}/reads/{{technology}}_{{accuracy}}/molephase/{region}.hap1.fq", region = define_phasing_regions()),
        expand(config["data_dir"] + "/{{coverage}}/reads/{{technology}}_{{accuracy}}/molephase/{region}.hap2.fq", region = define_phasing_regions()),
        expand(config["data_dir"] + "/{{coverage}}/reads/{{technology}}_{{accuracy}}/molephase/{region}.unphased.fq", region = define_phasing_regions()),
        fa_list = config["data_dir"] + "/{coverage}/global/fasta_list_{technology}_{accuracy}.tsv",
        binary = config["binary_dir"] + "create_consensus"
    output:
        consensus_dir = directory(config["data_dir"] + "/{coverage, [0-9]+}/global/consensus_seqs/{technology}_{accuracy}/{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}/{method}"),
        consensus_list = config["data_dir"] + "/{coverage, [0-9]+}/global/consensus_seqs_{method}_{technology}_{accuracy}_{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}.tsv",
        stats_out = config["data_dir"] + "/{coverage, [0-9]+}/global/consensus_stats_{method}_{technology}_{accuracy}_{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}.tsv"
    log: "logs/create_consensus_global_{method}_{coverage}_{technology}_{accuracy}_{p}.log"
    benchmark: "benchmarks/create_consensus_global_{method}_{coverage}_{technology}_{accuracy}_{p}.txt"
    threads: 16
    params:
        phasing_dir = lambda w: config["data_dir"] + "/" + w.coverage + "/reads/" + w.technology + "_" + w.accuracy + "/molephase",
        match = matchParam,
        mismatch = mismatchParam,
        gap = openParam,
        ext = extParam
    shell:
        "{input.binary} {input.fa_list} {output.consensus_dir} {output.consensus_list} {output.stats_out} {params.phasing_dir} {params.match} {params.mismatch} {params.gap} {params.ext} global {wildcards.method} &> {log}"


rule determine_consensus_times:
    input:
        expand(config["data_dir"] + "/{{coverage}}/reads/{{technology}}_{{accuracy}}/molephase/{region}.hap1.fq", region = define_phasing_regions()),
        expand(config["data_dir"] + "/{{coverage}}/reads/{{technology}}_{{accuracy}}/molephase/{region}.hap2.fq", region = define_phasing_regions()),
        expand(config["data_dir"] + "/{{coverage}}/reads/{{technology}}_{{accuracy}}/molephase/{region}.unphased.fq", region = define_phasing_regions()),
        fa_list = config["data_dir"] + "/{coverage}/fasta_list_{technology}_{accuracy}.tsv",
        binary = config["binary_dir"] + "determine_consensus_times"
    output:
        stats_out = config["data_dir"] + "/time_data/{coverage, [0-9]+}/consensus_stats_{method}_{technology}_{accuracy}_{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}.tsv"
    log: "logs/determine_consensus_times_{method}_{coverage}_{technology}_{accuracy}_{p}.log"
    benchmark: "benchmarks/determine_consensus_times_{method}_{coverage}_{technology}_{accuracy}_{p}.txt"
    threads: 1
    params:
        phasing_dir = lambda w: config["data_dir"] + "/" + w.coverage + "/reads/" + w.technology + "_" + w.accuracy + "/molephase",
        match = matchParam,
        mismatch = mismatchParam,
        gap = openParam,
        ext = extParam
    shell:
        "{input.binary} {input.fa_list} {output.stats_out} {params.phasing_dir} {params.match} {params.mismatch} {params.gap} {params.ext} {wildcards.method} &> {log}"

rule evaluate_consensus_sequences:
    input:
        consensus_list = config["data_dir"] + "/{coverage}/consensus_seqs_{method}_{technology}_{accuracy}_{p}.tsv",
        truth_dir = config["data_dir"] + "/allele_sequences",
        binary = config["binary_dir"] + "evaluate_consensus_seqs"
    output:
        config["results_dir"] + "/{technology}/{accuracy}/{coverage, [0-9]+}/{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}/evaluation_{method}.tsv"
    threads: 15
    log: "logs/eval_consensus_{method}_{coverage}_{technology}_{accuracy}_{p}.log"
    shell:
        '{input.binary} {input.consensus_list} {input.truth_dir} {output} &> {log}'

rule evaluate_consensus_sequences_semiglobal:
    input:
        consensus_list = config["data_dir"] + "/{coverage}/semiglobal/consensus_seqs_{method}_{technology}_{accuracy}_{p}.tsv",
        truth_dir = config["data_dir"] + "/allele_sequences",
        binary = config["binary_dir"] + "evaluate_consensus_seqs"
    output:
        config["results_dir"] + "/{technology}/{accuracy}/{coverage, [0-9]+}/{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}/evaluation_{method}_semiglobal.tsv"
    threads: 15
    log: "logs/eval_consensus_semiglobal_{method}_{coverage}_{technology}_{accuracy}_{p}.log"
    shell:
        '{input.binary} {input.consensus_list} {input.truth_dir} {output} &> {log}'

rule evaluate_consensus_sequences_spanning:
    input:
        consensus_list = config["data_dir"] + "/{coverage}/spanning/consensus_seqs_{method}_{technology}_{accuracy}_{p}.tsv",
        truth_dir = config["data_dir"] + "/allele_sequences",
        binary = config["binary_dir"] + "evaluate_consensus_seqs"
    output:
        config["results_dir"] + "/{technology}/{accuracy}/{coverage, [0-9]+}/{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}/evaluation_{method}_spanning.tsv"
    threads: 15
    log: "logs/eval_consensus_spanning_{method}_{coverage}_{technology}_{accuracy}_{p}.log"
    shell:
        '{input.binary} {input.consensus_list} {input.truth_dir} {output} &> {log}'

rule evaluate_consensus_sequences_global:
    input:
        consensus_list = config["data_dir"] + "/{coverage}/global/consensus_seqs_{method}_{technology}_{accuracy}_{p}.tsv",
        truth_dir = config["data_dir"] + "/allele_sequences",
        binary = config["binary_dir"] + "evaluate_consensus_seqs"
    output:
        config["results_dir"] + "/{technology}/{accuracy}/{coverage, [0-9]+}/{p, [0-9]+_-[0-9]+_-[0-9]+_-[0-9]+}/evaluation_{method}_global.tsv"
    threads: 15
    log: "logs/eval_consensus_global_{method}_{coverage}_{technology}_{accuracy}_{p}.log"
    shell:
        '{input.binary} {input.consensus_list} {input.truth_dir} {output} &> {log}'
