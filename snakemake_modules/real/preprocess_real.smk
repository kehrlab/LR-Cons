rule determine_hg002_regions:
    input:
        vcf = config["asm_dir"] + '/HG002_CHM13v2.0_v5.0q_stvar.vcf',
        sm_vcf = config["asm_dir"] + '/HG002_CHM13v2.0_v5.0q_smvar.vcf',
        ref = config["real_ref"]
    output:
        regions = config["asm_dir"] + '/variant_regions_{margin}.txt',
        sequences = directory(config["asm_dir"] + '/variant_sequences_{margin}')
    log:
        'logs/convert_HG002_vcf_to_regions_{margin}.log'
    resources:
        mem_mb = 16384
    conda: "../../envs/py_env.yaml"
    shell:
        '(python scripts/convert_vcf_to_regions_and_sim_giab_haplotypes.py {input.vcf} {input.sm_vcf} {input.ref} {wildcards.margin} {output.sequences} /dev/null > {output.regions}) 2> {log}'


rule mm2_idx_ont_real:
    input:
        fa = config["real_ref"]
    output:
        mmi = config["t2t_idx_ont"]
    conda: "../../envs/minimap_env.yaml"
    log: 'logs/index_t2t_ont.log'
    threads: 1
    resources:
        mem_mb = 30000
    shell:
        "minimap2 -x map-ont -d {output.mmi} {input.fa} &> {log}"

rule mm2_idx_pacbio_real:
    input:
        fa = config["real_ref"]
    output:
        mmi = config["t2t_idx_pb"]
    conda: "../../envs/minimap_env.yaml"
    log: 'logs/index_t2t_pacbio.log'
    threads: 1
    resources:
        mem_mb = 30000
    shell:
        "minimap2 -x map-pb -d {output.mmi} {input.fa} &> {log}"

rule align_pacbio_real:
    input:
        fq = config["realdata_dir"] + '/m84011_220902_175841_s1.hifi_reads.fastq',
        idx = config["t2t_idx_pb"]
    output:
        bam = temp(config["realdata_dir"] + '/t2t_remap/HG002_t2t_pacbio.bam'),
        bai = temp(config["realdata_dir"] + '/t2t_remap/HG002_t2t_pacbio.bam.bai')
    conda: "../../envs/minimap_env.yaml"
    log: 'logs/downloads/map_pb.log'
    threads:
        16
    params:
        p = "-a -Y --eqx -R '@RG\\tID:1\\tSM:HG002'"
    shell:
        "(minimap2 {params.p} -t {threads} {input.idx} {input.fq} | samtools sort -@ {threads} -o {output.bam} && samtools index {output.bam}) &> {log}"



rule realign_ont_real:
    input:
        bam = config["realdata_dir"] + '/ONT/{flowcell}/calls.sorted.bam',
        idx = config["t2t_idx_ont"]
    output:
        bam = temp(config["realdata_dir"] + '/ONT/{flowcell}/HG002_ont_realigned.bam')
    conda: "../../envs/minimap_env.yaml"
    log: 'logs/downloads/remap_ont_{flowcell}.log'
    threads:
        16
    params:
        p = lambda w:"-a -Y --eqx -R '@RG\\tID:" + w.flowcell + "\\tSM:HG002'"
    shell:
        "(samtools fasta {input.bam} | minimap2 {params.p} -t {threads} {input.idx} - | samtools sort -@ {threads} -o {output.bam}) &> {log}"

rule merge_ont_real:
    input:
        expand(config["realdata_dir"] + '/ONT/{flowcell}/HG002_ont_realigned.bam', flowcell = ['PAW70337', 'PAW71238'])
    output:
        bam = temp(config["realdata_dir"] + '/t2t_remap/HG002_t2t_ont.bam'),
        bai = temp(config["realdata_dir"] + '/t2t_remap/HG002_t2t_ont.bam.bai')
    conda: "../../envs/minimap_env.yaml"
    log:
        'logs/downloads/merge_ont.log'
    threads: 4
    shell:
        '(samtools merge -@ {threads} {output.bam} {input} && samtools index {output.bam}) &> {log}'


# --- subsample

rule subsample_real_bams:
    input:
        bam = config["realdata_dir"] + '/t2t_remap/HG002_t2t_{technology}.bam'
    output:
        bam = config["realdata_dir"] + '/t2t_remap/{coverage}/HG002_t2t_{technology}.bam',
        bai = config["realdata_dir"] + '/t2t_remap/{coverage}/HG002_t2t_{technology}.bam.bai'
    conda:
        "../../envs/minimap_env.yaml"
    log:
        "logs/subsample_{coverage}_{technology}.log"
    params:
        scaling = getScalingFactor # not the most efficient way, but keeps paths and rules consistent
    shell:
        '(samtools view{params.scaling} -b {input.bam} > {output.bam} && samtools index {output.bam}) &> {log}'


# a bit of QC
rule plot_coverage:
    input:
        bam = expand(config["realdata_dir"] + '/t2t_remap/{{coverage}}/HG002_t2t_{technology}.bam', technology = technologies),
        bai = expand(config["realdata_dir"] + '/t2t_remap/{{coverage}}/HG002_t2t_{technology}.bam.bai', technology = technologies)
    output:
        txt = config["realdata_dir"] + "/t2t_remap/{coverage}/coverages.tab",
        svg = config["realdata_dir"] + "/t2t_remap/{coverage}/coverages.svg"
    conda:
        "../../envs/lrqc.yaml"
    log: "logs/deeptools_plotCoverage_{coverage}.txt"
    priority: 1
    params: lambda w: '-n 100000 --plotTitle "Coverage Distribution ' + w.coverage + '" --ignoreDuplicates --plotHeight 10 --plotWidth 30 --plotFileFormat svg'
    shell:
        'plotCoverage -b {input.bam} {params} --outRawCounts {output.txt} --plotFile {output.svg} 2> {log}'


## --- sequence extraction and phasing
rule extract_real_region_sequences_real:
    input:
        bam = config["realdata_dir"] + "/t2t_remap/{coverage}/HG002_t2t_{technology}.bam",
        bai = config["realdata_dir"] + "/t2t_remap/{coverage}/HG002_t2t_{technology}.bam.bai",
        binary = config["binary_dir"] + "/extract_region_sequences",
        region_file = config["asm_dir"] + "/variant_regions_{margin}.txt"
    output:
        fa_list = config["realdata_dir"] + "/{mode}/fasta_list_{technology}_{coverage}_{margin}.tsv",
        dir = directory(config["realdata_dir"] + "/{mode}/region_sequences_{technology}_{coverage}_{margin}")
    log: "logs/extract_real_sequences_{technology}_{coverage}_{mode}_{margin}.log"
    shell:
        "{input.binary} {input.region_file} {output.dir} {output.fa_list} {input.bam} {wildcards.mode} &> {log}"

rule molephase_real:
    input:
        bam = config["realdata_dir"] + "/t2t_remap/{coverage}/HG002_t2t_{technology}.bam",
        bai = config["realdata_dir"] + "/t2t_remap/{coverage}/HG002_t2t_{technology}.bam.bai",
        ref = config["real_ref"],
        binary = "Molephase/molephase",
        region_file = config["asm_dir"] + "/variant_regions_1000.txt"
    output:
        directory(config["realdata_dir"] + "/reads/{technology}_{coverage}/molephase")
    log: "logs/molephase_real_{technology}_{coverage}.log"
    benchmark: "benchmarks/molephase/molephase_real_{technology}_{coverage}.txt"
    conda: "../../envs/molephase_env.yaml"
    threads: 8
    params:
        k = molephaseK
    shell:
        'python scripts/run_molephase_blocks.py {input.binary} {input.bam} {input.ref} {output} {params.k} {threads} {input.region_file} &> {log}'
