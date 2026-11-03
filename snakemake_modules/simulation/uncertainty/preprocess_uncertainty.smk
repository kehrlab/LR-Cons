rule extract_region_sequences_uncertainty:
    input:
        bam = config["simdata_dir"] + "/15/aln/{technology}/{technology}_reads_0.99.bam",
        bai = config["simdata_dir"] + "/15/aln/{technology}/{technology}_reads_0.99.bam.bai",
        binary = config["binary_dir"] + "/extract_region_sequences",
        region_file = config["simdata_dir"] + "/variant_regions_{margin}_{jitter}.tsv"
    output:
        fa_list = config["simdata_dir"] + "/uncertainty_{margin}_{jitter}/{mode}/fasta_list_{technology}.tsv",
        dir = directory(config["simdata_dir"] + "/uncertainty_{margin}_{jitter}/{mode}/region_sequences_{technology}.tsv")
    log: "logs/extract_sequences_{technology}_{mode}_{margin}_{jitter}.log"
    shell:
        "{input.binary} {input.region_file} {output.dir} {output.fa_list} {input.bam} {wildcards.mode} &> {log}"
