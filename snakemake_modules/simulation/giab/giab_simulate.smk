#
# Need to implement the python script called here
#
rule determine_hg002_regions_and_sim_giab_haplotypes:
    input:
        vcf = config["asm_dir"] + '/HG002_CHM13v2.0_v5.0q_stvar.vcf',
        sm_vcf = config["asm_dir"] + '/HG002_CHM13v2.0_v5.0q_smvar.vcf',
        ref = config["real_ref"]
    output:
        regions = config["simdata_dir"] + "/giab/variant_regions_500_0.tsv",
        sequences = directory(config["simdata_dir"] + "/giab/allele_sequences"),
        fa = config["simdata_dir"] + "/giab/haplotype_sequences.fa"
    log:
        'logs/convert_HG002_vcf_to_regions_and_sim_haplotypes_500_0.log'
    params:
        m = 500
    resources:
        mem_mb = 16384
    conda: "../../../envs/py_env.yaml"
    shell:
        '(python scripts/convert_vcf_to_regions_and_sim_giab_haplotypes.py {input.vcf} {input.sm_vcf} {input.ref} {params.m} {output.sequences} {output.fa} > {output.regions}) 2> {log}'

# simulate and align reads ...
#
rule simulate_ont_reads_giab:
    input:
        fa = config["simdata_dir"] + "/giab/haplotype_sequences.fa"
    output:
        fq = expand(config["simdata_dir"] + "/giab/ont/{{coverage}}/reads/ont_reads_{rg}.fq.gz", rg = ["0001", "0002"]),
        maf = temp(expand(config["simdata_dir"] + "/giab/ont/{{coverage}}/reads/ont_reads_{rg}.maf.gz", rg = ["0001", "0002"])),
        ref = temp(expand(config["simdata_dir"] + "/giab/ont/{{coverage}}/reads/ont_reads_{rg}.ref", rg = ["0001", "0002"]))
    conda:
        "../../../envs/pbsim_env.yaml"
    threads: 8
    resources:
        mem_mb = 4096
    log: "logs/giab_sim_pbsim_ont_{coverage}.log"
    params:
        model = config["pbsim_model_path"] + "/QSHMM-ONT-HQ.model",
        depth = lambda w: w.coverage,
        acc = lambda w: 0.99,
        prefix = lambda w: config["simdata_dir"] + "/giab/ont/" + w.coverage + "/reads/ont_reads"
    shell:
        "pbsim --strategy wgs --seed 34267861 --method qshmm --qshmm {params.model} --depth {params.depth} --accuracy-mean {params.acc} --genome {input.fa} --prefix {params.prefix} &> {log}"

rule merge_ont_fastq_giab:
    input:
        fq = expand(config["simdata_dir"] + "/giab/ont/{{coverage}}/reads/ont_reads_00{rg}.fq.gz", rg = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48"])
    output:
        fq = config["simdata_dir"] + "/giab/ont/{coverage}/reads/ont_reads.fq"
    log: 'logs/giab_sim_fastq_merge_{coverage}.log'
    resources:
        mem_mb = 8192
    shell:
        '(zcat {input.fq} > {output.fq}) &> {log}'

rule simulate_pacbio_reads_giab:
    input:
        fa = config["simdata_dir"] + "/giab/haplotype_sequences.fa"
    output:
        bam = expand(config["simdata_dir"] + "/giab/pacbio/{{coverage}}/reads/pacbio_reads_{{coverage}}_{rg}.bam", rg = ["0001", "0002"]),
        maf = temp(expand(config["simdata_dir"] + "/giab/pacbio/{{coverage}}/reads/pacbio_reads_{{coverage}}_{rg}.maf.gz", rg = ["0001", "0002"])),
        ref = temp(expand(config["simdata_dir"] + "/giab/pacbio/{{coverage}}/reads/pacbio_reads_{{coverage}}_{rg}.ref", rg = ["0001", "0002"]))
    log: "logs/giab_simp_bsim_pb_{coverage}.txt"
    benchmark: "benchmarks/giab_sim_pbsim_pb_{coverage}.txt"
    resources:
        mem_mb = 4096,
        runtime = 7200
    conda: "../../../envs/pbsim_env.yaml"
    params:
        model = config["pbsim_model_path"] + "/QSHMM-RSII.model",
        depth = lambda w: w.coverage,
        prefix = lambda w: config["simdata_dir"] + "/giab/pacbio/" + w.coverage + "/reads/pacbio_reads_" + w.coverage,
    shell:
        "(pbsim --strategy wgs --method qshmm --qshmm {params.model} --depth {params.depth} --genome {input.fa} --prefix {params.prefix} --pass-num 10) &> {log}"

rule pbccs_giab_sim:
    input:
        bam = config["simdata_dir"] + "/giab/pacbio/{coverage}/reads/pacbio_reads_{coverage}_00{rg}.bam"
    output:
        bam = config["simdata_dir"] + "/giab/pacbio/{coverage}/reads/pacbio_ccs_reads_{coverage}_00{rg}.bam"
    threads: 8
    resources:
        mem_mb = 16000,
        runtime = 7200
    conda: "../../../envs/pbsim_env.yaml"
    log: "logs/giab_sim_pbccs.{coverage}.{rg}.txt"
    benchmark: "benchmarks/giab_sim_pbccs.{coverage}.{rg}.txt"
    shell:
        "ccs --minPasses 3 --min-rq 0.99 -j {threads} {input.bam} {output.bam} &> {log}"
