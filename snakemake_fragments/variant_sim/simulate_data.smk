rule create_haplotype_sequences:
    input:
        ref = config["ref_file"],
        bed = config["exclude_file"],
        script = "scripts/create_haplotypes.py"
    output:
        fa = config["data_dir"] + "/chr1_haplotypes.fa",
        json = config["data_dir"] + "/chr1_var_info.json",
        truth_dir = directory(config["data_dir"] + "/allele_sequences")
    conda: "../../envs/py_env.yaml"
    shell:
        'python {input.script} {input.ref} {input.bed} chr1 {output.json} {output.fa} {output.truth_dir}'

rule simulate_ont_reads:
    input:
        fa = config["data_dir"] + "/chr1_haplotypes.fa"
    output:
        fq = expand(config["data_dir"] + "/{{coverage, [0-9]+}}/reads/ont/ont_reads_{{accuracy}}_{rg}.fq.gz", rg = ["0001", "0002"]),
        maf = temp(expand(config["data_dir"] + "/{{coverage, [0-9]+}}/reads/ont/ont_reads_{{accuracy}}_{rg}.maf.gz", rg = ["0001", "0002"])),
        ref = temp(expand(config["data_dir"] + "/{{coverage, [0-9]+}}/reads/ont/ont_reads_{{accuracy}}_{rg}.ref", rg = ["0001", "0002"]))
    conda:
        "../../envs/pbsim_env.yaml"
    threads: 8
    log: "logs/pbsim_ont_{coverage}_{accuracy}.log"
    params:
        model = config["pbsim_model_path"] + "/QSHMM-ONT-HQ.model",
        depth = lambda w: w.coverage,
        acc = lambda w: w.accuracy,
        prefix = lambda w: config["data_dir"] + "/" + w.coverage + "/reads/ont/ont_reads_" + w.accuracy
    shell:
        "pbsim --strategy wgs --seed 34267861 --method qshmm --qshmm {params.model} --depth {params.depth} --accuracy-mean {params.acc} --genome {input.fa} --prefix {params.prefix} &> {log}"

rule merge_ont_fastq:
    input:
        fq1 = config["data_dir"] + "/{coverage}/reads/ont/ont_reads_{accuracy}_0001.fq.gz",
        fq2 = config["data_dir"] + "/{coverage}/reads/ont/ont_reads_{accuracy}_0002.fq.gz"
    output:
        fq = config["data_dir"] + "/{coverage, [0-9]+}/reads/ont/ont_reads_{accuracy}.fq"
    log: 'logs/fastq_merge_{coverage}_{accuracy}.log'
    shell:
        '(zcat {input.fq1} > {output.fq} && zcat {input.fq2} >> {output.fq}) &> {log}'

rule simulate_pacbio_reads:
    input:
        fa = config["data_dir"] + "/chr1_haplotypes.fa"
    output:
        bam = expand(config["data_dir"] + "/{{coverage, [0-9]+}}/reads/pacbio/pacbio_reads_{{coverage}}_{rg}.bam", rg = ["0001", "0002"]),
        maf = temp(expand(config["data_dir"] + "/{{coverage, [0-9]+}}/reads/pacbio/pacbio_reads_{{coverage}}_{rg}.maf.gz", rg = ["0001", "0002"])),
        ref = temp(expand(config["data_dir"] + "/{{coverage, [0-9]+}}/reads/pacbio/pacbio_reads_{{coverage}}_{rg}.ref", rg = ["0001", "0002"]))
    log: "logs/pbsim_pb_{coverage}.txt"
    benchmark: "benchmarks/pbsim_pb_{coverage}.txt"
    conda: "../../envs/pbsim_env.yaml"
    params:
        model = config["pbsim_model_path"] + "/QSHMM-RSII.model",
        depth = lambda w: w.coverage,
        prefix = lambda w: config["data_dir"] + "/" + w.coverage + "/reads/pacbio/pacbio_reads_" + w.coverage,
    shell:
        "(pbsim --strategy wgs --method qshmm --qshmm {params.model} --depth {params.depth} --genome {input.fa} --prefix {params.prefix} --pass-num 10) &> {log}"

rule pbccs:
    input:
        bam = config["data_dir"] + "/{coverage, [0-9]+}/reads/pacbio/pacbio_reads_{coverage}_{rg}.bam"
    output:
        bam = config["data_dir"] + "/{coverage, [0-9]+}/reads/pacbio/pacbio_ccs_reads_{coverage}_{rg}.bam"
    threads: 8
    conda: "../../envs/pbsim_env.yaml"
    log: "logs/pbccs.{coverage}.{rg}.txt"
    benchmark: "benchmarks/pbccs.{coverage}.{rg}.txt"
    shell:
        "ccs --minPasses 3 --min-rq 0.99 -j {threads} {input.bam} {output.bam} &> {log}"


rule define_regions:
    input:
        json = config["data_dir"] + "/chr1_var_info.json",
        script = "scripts/define_regions.py"
    output:
        config["data_dir"] + "/variant_regions.tsv"
    conda: "../../envs/py_env.yaml"
    log: "logs/define_regions.log"
    params:
        margin = 500
    shell:
        "python {input.script} {input.json} {params.margin} {output} &> {log}"
