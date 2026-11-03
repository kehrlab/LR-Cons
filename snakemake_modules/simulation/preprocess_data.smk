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
        fq = config["simdata_dir"] + "/{coverage}/reads/ont/ont_reads_{accuracy}.fq"
    output:
        bam = config["simdata_dir"] + "/{coverage, [0-9]+}/aln/ont/ont_reads_{accuracy}.bam"
    conda: "../../envs/minimap_env.yaml"
    log: "logs/mm2_ont_{coverage}_{accuracy}.log"
    threads: 8
    resources:
        mem_mb = 24000
    shell:
        "(minimap2 -a -Y -t {threads} --rmq=yes -R @RG\\\\tID:0\\\\tSM:Sample_0 {input.idx} {input.fq} | samtools sort -@ {threads} -O bam -o {output.bam}) &> {log}"#@RG\\\\tID:{wildcards.rg}\\\\tSM:Sample_0

rule mm2_map_hifi:
    input:
        bam = config["simdata_dir"] + "/{coverage}/reads/pacbio/pacbio_ccs_reads_{coverage}_{rg}.bam",
        idx = config["mm_idx_pb"]
    output:
        bam = temp(config["simdata_dir"] + "/{coverage}/aln/pacbio/pacbio_ccs_reads_{rg}.bam")
    log: "logs/minimap2_ccs_{coverage}_{rg}.txt"
    benchmark: "benchmarks/minimap2_ccs_{coverage}_{rg}.txt"
    threads: 8
    resources:
        mem_mb = 24000
    conda: "../../envs/minimap_env.yaml"
    shell:
        "(samtools fastq -n -t {input.bam} | minimap2 -a -Y -t {threads} --rmq=yes -R @RG\\\\tID:{wildcards.rg}\\\\tSM:Sample_0 {input.idx} - | samtools sort -@ {threads} -O bam -o {output.bam}) &> {log}"


rule merge_pacbio:
    input:
        bam = expand(config["simdata_dir"] + "/{{coverage}}/aln/pacbio/pacbio_ccs_reads_{rg}.bam", rg = ["0001", "0002"])
    output:
        bam = config["simdata_dir"] + "/{coverage}/aln/pacbio/pacbio_reads_{accuracy}.bam"
    conda: "../../envs/minimap_env.yaml"
    threads: 8
    resources:
        mem_mb = 4096
    log: "logs/merge_bp_{coverage}_{accuracy}.log"
    shell:
        "samtools merge -@ {threads} {output.bam} {input.bam} &> {log}"

rule index:
    input:
        bam = config["simdata_dir"] + "/{coverage}/aln/{technology}/{technology}_reads_{accuracy}.bam"
    output:
        idx = config["simdata_dir"] + "/{coverage, [0-9]+}/aln/{technology}/{technology}_reads_{accuracy}.bam.bai"
    conda: "../../envs/minimap_env.yaml"
    threads: 8
    log: "logs/index_{technology}_{coverage}_{accuracy}.log"
    shell:
        "samtools index {input.bam} > {log}"

rule molephase:
    input:
        bam = config["simdata_dir"] + "/{coverage}/aln/{technology}/{technology}_reads_{accuracy}.bam",
        bai = config["simdata_dir"] + "/{coverage}/aln/{technology}/{technology}_reads_{accuracy}.bam.bai",
        ref = config["ref_file"],
        binary = "Molephase/molephase"
    output:
        directory(config["simdata_dir"] + "/{coverage, [0-9]+}/reads/{technology}_{accuracy}/molephase")
    log: "logs/molephase_{coverage}_{technology}_{accuracy}.log"
    benchmark: "benchmarks/molephase/molephase_{technology}_{accuracy}_{coverage}.txt"
    conda: "../../envs/molephase_env.yaml"
    threads: 1
    params:
        k = molephaseK
    shell:
        'python scripts/run_molephase_blocks.py {input.binary} {input.bam} {input.ref} {output} {params.k} {threads} &> {log}'
