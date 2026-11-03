# --- HG002 PacBio download
rule download_pacbio_hg002_old:
    output:
        fq = temp(config["realdata_dir"] + '/PB_old/{sample}')
    log:
        'logs/downloads/hg002_pb_old_{sample}.log'
    shell:
        'wget -O {output.fq} ' + config["pb_old_url"] + '{wildcards.sample} 2> {log}'

#rule download_pacbio_hg002_old:
#    output:
#        fq2 = temp(config["realdata_dir"] + '/PB_old/m64012_190921_234837.Q20.fastq'),
#        fq1 = temp(config["realdata_dir"] + '/PB_old/m64012_190920_173625.Q20.fastq'),
#        fq3 = temp(config["realdata_dir"] + '/PB_old/m64011_190830_220126.Q20.fastq'),
#        fq4 = temp(config["realdata_dir"] + '/PB_old/m64011_190901_095311.Q20.fastq')
#    log:
#        'logs/downloads/hg002_pb_old.log'
#    shell:
#        'wget -O {output.fq1} ' + config["pb_old_url_1"] + ' && wget -O {output.fq2} ' + config["pb_old_url_2"] + ' && wget -O {output.fq3} ' + config["pb_old_url_3"] + ' &&  wget -O {output.fq4} ' + config["pb_old_url_4"] + ' &'


# --- HG002 ONT download
rule download_ont_hg002_old:
    output:
        temp(config["realdata_dir"] + '/ONT_old/basecalls.fastq.gz')
    conda:
        '../../envs/aws_env.yaml'
    log:
        'logs/downloads/hg002_ont_old.log'
    shell:
        'aws s3 --no-sign-request cp ' + config["ont_aws_old"] + ' {output} &> {log}'
