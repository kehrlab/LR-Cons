rule mm2_map_ont_giab_sim:
    input:
        idx = config["t2t_idx_ont"],
        fq = config["simdata_dir"] + "/giab/ont/{coverage}/reads/ont_reads.fq"
    output:
        bam = config["simdata_dir"] + "/giab/ont/{coverage}/aln/ont_reads.bam"
    conda: "../../../envs/minimap_env.yaml"
    log: "logs/mm2_ont_{coverage}.log"
    threads: 12
    resources:
        mem_mb = 48000,
        cpus_per_task = 12
    shell:
        "(minimap2 -a -Y -t {threads} --rmq=yes -R @RG\\\\tID:0\\\\tSM:Sample_0 {input.idx} {input.fq} | samtools sort -@ {threads} -O bam -o {output.bam}) &> {log}"#@RG\\\\tID:{wildcards.rg}\\\\tSM:Sample_0

rule mm2_map_hifi_giab_sim:
    input:
        bam = config["simdata_dir"] + "/giab/pacbio/{coverage}/reads/pacbio_ccs_reads_{coverage}_00{rg}.bam",
        idx = config["t2t_idx_pb"]
    output:
        bam = config["simdata_dir"] + "/giab/pacbio/{coverage}/aln/pacbio_ccs_reads_00{rg}.bam"
    log: "logs/giab_sim_minimap2_ccs_{coverage}_{rg}.txt"
    benchmark: "benchmarks/giab_sim_minimap2_ccs_{coverage}_{rg}.txt"
    threads: 8
    resources:
        mem_mb = 24000
    conda: "../../../envs/minimap_env.yaml"
    shell:
        "(samtools fastq -n -t {input.bam} | minimap2 -a -Y -t {threads} --rmq=yes -R @RG\\\\tID:{wildcards.rg}\\\\tSM:Sample_0 {input.idx} - | samtools sort -@ {threads} -O bam -o {output.bam}) &> {log}"


rule merge_pacbio_giab_sim:
    input:
        bam = expand(config["simdata_dir"] + "/giab/pacbio/{{coverage}}/aln/pacbio_ccs_reads_00{rg}.bam", rg = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48"])
    output:
        bam = config["simdata_dir"] + "/giab/pacbio/{coverage}/aln/pacbio_reads.bam"
    conda: "../../../envs/minimap_env.yaml"
    threads: 8
    resources:
        mem_mb = 4096
    log: "logs/giab_sim_merge_bp_{coverage}.log"
    shell:
        "samtools merge -@ {threads} {output.bam} {input.bam} &> {log}"

rule index_giab_sim:
    input:
        bam = config["simdata_dir"] + "/giab/{technology}/{coverage}/aln/{technology}_reads.bam"
    output:
        idx = config["simdata_dir"] + "/giab/{technology}/{coverage}/aln/{technology}_reads.bam.bai"
    conda: "../../../envs/minimap_env.yaml"
    threads: 8
    log: "logs/giab_sim_index_{technology}_{coverage}.log"
    shell:
        "samtools index {input.bam} > {log}"


rule extract_region_sequences_giab_sim:
    input:
        bam = config["simdata_dir"] + "/giab/{technology}/{coverage}/aln/{technology}_reads.bam",
        bai = config["simdata_dir"] + "/giab/{technology}/{coverage}/aln/{technology}_reads.bam.bai",
        binary = config["binary_dir"] + "/extract_region_sequences",
        region_file = config["simdata_dir"] + "/giab/variant_regions_500_0.tsv"
    output:
        fa_list = config["simdata_dir"] + "/giab/{technology}/{coverage}/fasta_list_{technology}_{mode}.tsv",
        dir = directory(config["simdata_dir"] + "/giab/{technology}/{coverage}/region_sequences/{mode}")
    log: "logs/giab_sim_extract_sequences_{coverage}_{technology}_{mode}.log"
    shell:
        "{input.binary} {input.region_file} {output.dir} {output.fa_list} {input.bam} {wildcards.mode} &> {log}"


rule molephase_giab_sim:
    input:
        bam = config["simdata_dir"] + "/giab/{technology}/{coverage}/aln/{technology}_reads.bam",
        bai = config["simdata_dir"] + "/giab/{technology}/{coverage}/aln/{technology}_reads.bam.bai",
        ref = config["real_ref"],
        binary = "Molephase/molephase",
        region_file = config["simdata_dir"] + "/giab/variant_regions_500_0.tsv"
    output:
        directory(config["simdata_dir"] + "/giab/{technology}/{coverage}/reads/{technology}_0.99/molephase")
    log: "logs/giab_sim_molephase_{coverage}_{technology}.log"
    benchmark: "benchmarks/molephase/giab_sim_molephase_{technology}_{coverage}.txt"
    conda: "../../../envs/molephase_env.yaml"
    threads: 8
    priority: 5
    params:
        k = molephaseK
    resources:
        mem_mb = 8192,
        cpus_per_task = 8
    shell:
        'python scripts/run_molephase_blocks.py {input.binary} {input.bam} {input.ref} {output} {params.k} {threads} {input.region_file} &> {log}'
