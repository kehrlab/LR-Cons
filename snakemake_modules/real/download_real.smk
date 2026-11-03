# --- HG002 assembly and VCF download (w.r.t. T2T-CHM13)

rule download_hg002_assembly:
    output:
        fa = config["asm_dir"] + '/hg002v1.1.fasta'
    log:
        'logs/downloads/hg002_assembly.log'
    shell:
        '(wget -O {output.fa}.gz ' + config["asm_url"] + " && gunzip {output.fa}.gz) &> {log}"

rule download_giab_vcf:
    output:
        vcf = config["asm_dir"] + '/HG002_CHM13v2.0_v5.0q_stvar.vcf'
    log:
        'logs/downloads/vcf.log'
    shell:
        '(wget -O {output.vcf}.gz ' + config["vcf_url"] + " && gunzip {output.vcf}.gz) &> {log}"

rule download_giab_small_vcf:
    output:
        vcf = config["asm_dir"] + '/HG002_CHM13v2.0_v5.0q_smvar.vcf'
    log:
        'logs/downloads/vcf_sm.log'
    shell:
        '(wget -O {output.vcf}.gz ' + config["sm_vcf_url"] + " && gunzip {output.vcf}.gz) &> {log}"

# --- HG002 PacBio download
rule download_pacbio_hg002:
    output:
        fq = temp(config["realdata_dir"] + '/m84011_220902_175841_s1.hifi_reads.fastq')
    log:
        'logs/downloads/hg002_pb.log'
    shell:
        '(wget -O {output.fq}.gz ' + config["pb_url"] + ' && gunzip {output.fq}.gz) &> {log}'

# --- HG002 ONT download
rule download_ont_hg002:
    output:
        bam1 = temp(config["realdata_dir"] + '/ONT/PAW70337/calls.sorted.bam'),
        bam2 = temp(config["realdata_dir"] + '/ONT/PAW71238/calls.sorted.bam')
    conda:
        '../../envs/aws_env.yaml'
    log:
        'logs/downloads/hg002_ont.log'
    shell:
        'aws s3 sync --no-sign-request ' + config["ont_aws"] + ' ' + config["realdata_dir"] + '/ONT/ &> {log}'

rule download_hg38:
    output:
        fa = 'resources/hg38.fna'
    log:
        'logs/downloads/hg38.log'
    shell:
        '(wget -O {output.fa}.gz ' + config['hg38_url'] + ' && gunzip {output.fa}.gz) &> {log}'
