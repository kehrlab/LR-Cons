#rule determine_hg002_regions_old:
#    input:
#        vcf = config["asm_dir"] + '/HG002_CHM13v2.0_v5.0q_stvar.vcf',
#        sm_vcf = config["asm_dir"] + '/HG002_CHM13v2.0_v5.0q_smvar.vcf',
#        ref = config["real_ref"]
#    output:
#        regions = config["asm_dir"] + '/variant_regions_{margin}.txt',
#        sequences = directory(config["asm_dir"] + '/variant_sequences_{margin}')
#    log:
#        'logs/convert_HG002_vcf_to_regions_{margin}.log'
#    resources:
#        mem_mb = 16384
#    conda: "../../envs/py_env.yaml"
#    shell:
#        '(python scripts/convert_vcf_to_regions_and_sim_giab_haplotypes.py {input.vcf} {input.sm_vcf} {input.ref} {wildcards.margin} {output.sequences} /dev/null > {output.regions}) 2> {log}'


rule merge_pacbio_real_old:
    input:
        fq1 = config["realdata_dir"] + '/PB_old/15kb/m64012_190920_173625.Q20.fastq',
        fq2 = config["realdata_dir"] + '/PB_old/15kb/m64012_190921_234837.Q20.fastq',
        fq3 = config["realdata_dir"] + '/PB_old/20kb/m64011_190830_220126.Q20.fastq',
        fq4 = config["realdata_dir"] + '/PB_old/20kb/m64011_190901_095311.Q20.fastq'
    output:
        fq = temp(config["realdata_dir"] + '/PB_old/basecalls.fastq')
    log: 'logs/merge_pb_real_old.log'
    threads: 1
    resources:
        mem_mb=16384
    shell: 
        "cat {input.fq1} {input.fq2} {input.fq3} {input.fq4} > {output.fq} 2> {log}"

rule align_pacbio_real_old:
    input:
        fq = config["realdata_dir"] + '/PB_old/basecalls.fastq',
        idx = config["t2t_idx_pb"]
    output:
        bam = config["realdata_dir"] + '/t2t_remap/HG002_t2t_pacbio_old.bam',
        bai = config["realdata_dir"] + '/t2t_remap/HG002_t2t_pacbio_old.bam.bai'
    conda: "../../envs/minimap_env.yaml"
    log: 'logs/downloads/map_pb_old.log'
    resources:
        mem_mb = 48000,
        cpus_per_task = 16
    threads:
        16
    params:
        p = "-a -Y --eqx -R '@RG\\tID:1\\tSM:HG002'"
    shell:
        "(minimap2 {params.p} -t {threads} {input.idx} {input.fq} | samtools sort -@ {threads} -o {output.bam} && samtools index {output.bam}) &> {log}"



rule align_ont_real_old:
    input:
        fq = config["realdata_dir"] + '/ONT_old/basecalls.fastq.gz',
        idx = config["t2t_idx_ont"]
    output:
        bam = config["realdata_dir"] + '/t2t_remap/HG002_t2t_ont_old.bam'
    conda: "../../envs/minimap_env.yaml"
    log: 'logs/downloads/remap_ont_old.log'
    threads:
        16
    resources:
        mem_mb = 48000,
        cpus_per_task = 16
    params:
        p = lambda w:"-a -Y --eqx -R '@RG\\tID:1\\tSM:HG002'"
    shell:
        "ulimit -n 2048; (zcat {input.fq} | minimap2 {params.p} -t {threads} {input.idx} - | samtools sort -@ {threads} -o {output.bam}) &> {log}"

#rule merge_ont_real:
#    input:
#        config["realdata_dir"] + '/ONT_old/HG002_ont_old_realigned.bam'
#    output:
#        bam = temp(config["realdata_dir"] + '/t2t_remap/HG002_t2t_ont_old.bam'),
#        bai = temp(config["realdata_dir"] + '/t2t_remap/HG002_t2t_ont_old.bam.bai')
#    conda: "../../envs/minimap_env.yaml"
#    log:
#        'logs/downloads/merge_ont_old.log'
#    threads: 4
#    shell:
#        '(samtools merge -@ {threads} {output.bam} {input} && samtools index {output.bam}) &> {log}'


# --- subsample

rule subsample_real_bams_old:
    input:
        bam = config["realdata_dir"] + '/t2t_remap/HG002_t2t_{technology}_old.bam'
    output:
        bam = config["realdata_dir"] + '/t2t_remap/{coverage}/HG002_t2t_{technology}_old.bam',
        bai = config["realdata_dir"] + '/t2t_remap/{coverage}/HG002_t2t_{technology}_old.bam.bai'
    conda:
        "../../envs/minimap_env.yaml"
    log:
        "logs/subsample_{coverage}_{technology}_old.log"
    params:
        scaling = getScalingFactorOld # not the most efficient way, but keeps paths and rules consistent
    shell:
        '(samtools view{params.scaling} -b {input.bam} > {output.bam} && samtools index {output.bam}) &> {log}'


# a bit of QC
rule plot_coverage_old:
    input:
        bam = expand(config["realdata_dir"] + '/t2t_remap/{{coverage}}/HG002_t2t_{technology}_old.bam', technology = technologies),
        bai = expand(config["realdata_dir"] + '/t2t_remap/{{coverage}}/HG002_t2t_{technology}_old.bam.bai', technology = technologies)
    output:
        txt = config["realdata_dir"] + "/t2t_remap/{coverage}/coverages_old.tab",
        svg = config["realdata_dir"] + "/t2t_remap/{coverage}/coverages_old.svg"
    conda:
        "../../envs/lrqc.yaml"
    log: "logs/deeptools_plotCoverage_{coverage}_old.txt"
    priority: 1
    params: lambda w: '-n 100000 --plotTitle "Coverage Distribution ' + w.coverage + '" --ignoreDuplicates --plotHeight 10 --plotWidth 30 --plotFileFormat svg'
    shell:
        'plotCoverage -b {input.bam} {params} --outRawCounts {output.txt} --plotFile {output.svg} 2> {log}'


## --- sequence extraction and phasing
rule extract_real_region_sequences_real_old:
    input:
        bam = config["realdata_dir"] + "/t2t_remap/{coverage}/HG002_t2t_{technology}_old.bam",
        bai = config["realdata_dir"] + "/t2t_remap/{coverage}/HG002_t2t_{technology}_old.bam.bai",
        binary = config["binary_dir"] + "/extract_region_sequences",
        region_file = config["asm_dir"] + "/variant_regions_{margin}.txt"
    output:
        fa_list = config["realdata_dir"] + "/old/{mode}/fasta_list_{technology}_{coverage}_{margin}.tsv",
        dir = directory(config["realdata_dir"] + "/old/{mode}/region_sequences_{technology}_{coverage}_{margin}")
    log: "logs/extract_real_sequences_{technology}_{coverage}_{mode}_{margin}_old.log"
    shell:
        "{input.binary} {input.region_file} {output.dir} {output.fa_list} {input.bam} {wildcards.mode} &> {log}"

rule molephase_real_old:
    input:
        bam = config["realdata_dir"] + "/t2t_remap/{coverage}/HG002_t2t_{technology}_old.bam",
        bai = config["realdata_dir"] + "/t2t_remap/{coverage}/HG002_t2t_{technology}_old.bam.bai",
        ref = config["real_ref"],
        binary = "Molephase/molephase",
        region_file = config["asm_dir"] + "/variant_regions_1000.txt"
    output:
        directory(config["realdata_dir"] + "/old/reads/{technology}_{coverage}/molephase")
    log: "logs/molephase_real_{technology}_{coverage}_old.log"
    benchmark: "benchmarks/molephase/molephase_real_{technology}_{coverage}_old.txt"
    conda: "../../envs/molephase_env.yaml"
    threads: 8
    resources:
        mem_mb=8192,
        cpus_per_task=8
    params:
        k = molephaseK
    shell:
        'python scripts/run_molephase_blocks.py {input.binary} {input.bam} {input.ref} {output} {params.k} {threads} {input.region_file} &> {log}'
