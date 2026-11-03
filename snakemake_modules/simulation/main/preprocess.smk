rule extract_region_sequences:
    input:
        bam = config["simdata_dir"] + "/{coverage}/aln/{technology}/{technology}_reads_{accuracy}.bam",
        bai = config["simdata_dir"] + "/{coverage}/aln/{technology}/{technology}_reads_{accuracy}.bam.bai",
        binary = config["binary_dir"] + "/extract_region_sequences",
        region_file = config["simdata_dir"] + "/variant_regions_500_0.tsv"
    output:
        fa_list = config["simdata_dir"] + "/{coverage, [0-9]+}/{mode}/fasta_list_{technology}_{accuracy}.tsv",
        dir = directory(config["simdata_dir"] + "/{coverage, [0-9]+}/{mode}/region_sequences_{technology}_{accuracy}.tsv")
    log: "logs/extract_sequences_{coverage}_{technology}_{accuracy}_{mode}.log"
    shell:
        "{input.binary} {input.region_file} {output.dir} {output.fa_list} {input.bam} {wildcards.mode} &> {log}"
