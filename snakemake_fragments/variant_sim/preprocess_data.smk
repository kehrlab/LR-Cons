rule mm2_idx_ont:
    input:
        ref = config["ref_file"]
    output:
        idx = config["mm_idx_ont"]
    threads: 8
    conda: "../../envs/minimap_env.yaml"
    shell:
        "minimap2 -t {threads} -x map-ont -d {output.idx} {input.ref}"

rule mm2_idx_pacbio:
    input:
        ref = config["ref_file"]
    output:
        idx = config["mm_idx_pb"]
    threads: 8
    conda: "../../envs/minimap_env.yaml"
    shell:
        "minimap2 -t {threads} -x map-hifi -d {output.idx} {input.ref}"

rule mm2_map_ont:
    input:
        idx = config["mm_idx_ont"],
        fq = config["data_dir"] + "/{coverage}/reads/ont/ont_reads_{accuracy}.fq"
    output:
        bam = config["data_dir"] + "/{coverage, [0-9]+}/aln/ont/ont_reads_{accuracy}.bam"
    conda: "../../envs/minimap_env.yaml"
    log: "logs/mm2_ont_{coverage}_{accuracy}.log"
    threads: 16
    shell:
        "(minimap2 -a -Y -t {threads} --rmq=yes -R @RG\\\\tID:0\\\\tSM:Sample_0 {input.idx} {input.fq} | samtools sort -@ {threads} -O bam -o {output.bam}) &> {log}"#@RG\\\\tID:{wildcards.rg}\\\\tSM:Sample_0

rule mm2_map_hifi:
    input:
        bam = config["data_dir"] + "/{coverage}/reads/pacbio/pacbio_ccs_reads_{coverage}_{rg}.bam",
        idx = config["mm_idx_pb"]
    output:
        bam = temp(config["data_dir"] + "/{coverage}/aln/pacbio/pacbio_ccs_reads_{rg}.bam")
    log: "logs/minimap2_ccs_{coverage}_{rg}.txt"
    benchmark: "benchmarks/minimap2_ccs_{coverage}_{rg}.txt"
    threads: 16
    conda: "../../envs/minimap_env.yaml"
    shell:
        "(samtools fastq -n -t {input.bam} | minimap2 -a -Y -t {threads} --rmq=yes -R @RG\\\\tID:{wildcards.rg}\\\\tSM:Sample_0 {input.idx} - | samtools sort -@ {threads} -O bam -o {output.bam}) &> {log}"


rule merge_pacbio:
    input:
        bam = expand(config["data_dir"] + "/{{coverage}}/aln/pacbio/pacbio_ccs_reads_{rg}.bam", rg = ["0001", "0002"])
    output:
        bam = config["data_dir"] + "/{coverage}/aln/pacbio/pacbio_reads_{accuracy}.bam"
    conda: "../../envs/minimap_env.yaml"
    threads: 8
    log: "logs/merge_bp_{coverage}_{accuracy}.log"
    shell:
        "samtools merge -@ {threads} {output.bam} {input.bam} &> {log}"

rule index:
    input:
        bam = config["data_dir"] + "/{coverage}/aln/{technology}/{technology}_reads_{accuracy}.bam"
    output:
        idx = config["data_dir"] + "/{coverage, [0-9]+}/aln/{technology}/{technology}_reads_{accuracy}.bam.bai"
    conda: "../../envs/minimap_env.yaml"
    threads: 8
    log: "logs/index_{technology}_{coverage}_{accuracy}.log"
    shell:
        "samtools index {input.bam} > {log}"

rule extract_region_sequences:
    input:
        bam = config["data_dir"] + "/{coverage}/aln/{technology}/{technology}_reads_{accuracy}.bam",
        bai = config["data_dir"] + "/{coverage}/aln/{technology}/{technology}_reads_{accuracy}.bam.bai",
        binary = config["binary_dir"] + "extract_region_sequences",
        region_file = config["data_dir"] + "/variant_regions.tsv"
    output:
        fa_list = config["data_dir"] + "/{coverage, [0-9]+}/fasta_list_{technology}_{accuracy}.tsv",
        dir = directory(config["data_dir"] + "/{coverage, [0-9]+}/region_sequences_{technology}_{accuracy}")
    log: "logs/extract_sequences_{coverage}_{technology}_{accuracy}.log"
    shell:
        "{input.binary} {input.region_file} {output.dir} {output.fa_list} {input.bam} PARTIAL &> {log}"

rule extract_region_sequences_global:
    input:
        bam = config["data_dir"] + "/{coverage}/aln/{technology}/{technology}_reads_{accuracy}.bam",
        bai = config["data_dir"] + "/{coverage}/aln/{technology}/{technology}_reads_{accuracy}.bam.bai",
        binary = config["binary_dir"] + "extract_region_sequences",
        region_file = config["data_dir"] + "/variant_regions.tsv"
    output:
        fa_list = config["data_dir"] + "/{coverage, [0-9]+}/global/fasta_list_{technology}_{accuracy}.tsv",
        dir = directory(config["data_dir"] + "/{coverage, [0-9]+}/global/region_sequences_{technology}_{accuracy}")
    log: "logs/extract_sequences_global_{coverage}_{technology}_{accuracy}.log"
    shell:
        "{input.binary} {input.region_file} {output.dir} {output.fa_list} {input.bam} COMPLETE &> {log}"

rule molephase:
    input:
        bam = config["data_dir"] + "/{coverage}/aln/{technology}/{technology}_reads_{accuracy}.bam",
        bai = config["data_dir"] + "/{coverage}/aln/{technology}/{technology}_reads_{accuracy}.bam.bai",
        ref = config["ref_file"],
        binary = "Molephase/molephase"
    output:
        hap1 = config["data_dir"] + "/{coverage, [0-9]+}/reads/{technology}_{accuracy}/molephase/{region}.hap1.fq",
        hap2 = config["data_dir"] + "/{coverage, [0-9]+}/reads/{technology}_{accuracy}/molephase/{region}.hap2.fq",
        unphased = config["data_dir"] + "/{coverage, [0-9]+}/reads/{technology}_{accuracy}/molephase/{region}.unphased.fq",
        all = temp(config["data_dir"] + "/{coverage, [0-9]+}/reads/{technology}_{accuracy}/molephase/{region}_all.fq")
    log: "logs/molephase_{region}_{coverage}_{technology}_{accuracy}.log"
    benchmark: "benchmarks/molephase/molephase_{region}_{coverage}_{technology}_{accuracy}.txt"
    conda: "../../envs/minimap_env.yaml"
    params:
        start = startParam,
        end = endParam,
        chr = chrParam,
        prefix = lambda w: config["data_dir"] + "/" + w.coverage + "/reads/" + w.technology + "_" + w.accuracy + "/molephase/" + w.region,
        k = molephaseK
    shell:
        '(samtools view -h {input.bam} {params.chr}:{params.start}-{params.end} | samtools fastq - > {output.all}; ./{input.binary} -1 {output.all} -o {params.prefix} -r {input.ref} -c {params.chr} -s {params.start} -e {params.end} -k {params.k}) &> {log}'
