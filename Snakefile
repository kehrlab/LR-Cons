configfile: "config.yaml"

# --- helper functions for parameter extraction
def chrParam(wildcards):
    return wildcards.region.split("_")[0]

def startParam(wildcards):
    return wildcards.region.split("_")[1]

def endParam(wildcards):
    return wildcards.region.split("_")[2]

def matchParam(wildcards):
    return wildcards.p.split("_")[0]

def mismatchParam(wildcards):
    return wildcards.p.split("_")[1]

def openParam(wildcards):
    return wildcards.p.split("_")[2]

def extParam(wildcards):
    return wildcards.p.split("_")[3]

def molephaseK(wildcards):
    if (wildcards.technology == "ont"):
        return 55
    elif (wildcards.technology == "pacbio"):
        return 77

# --- helper functions for
def getFullCov(wildcards):
    if wildcards.technology == "ont":
        return 94.10
    elif wildcards.technology == "pacbio":
        return 32.51

def getScalingFactor(wildcards):
    seed = "862026"
    if wildcards.technology == "ont":
        if wildcards.coverage == 94.1:
            return ""
        else:
            return " -s " + seed + "." + str(float(wildcards.coverage) / 94.1)[2:]
    elif wildcards.technology == "pacbio":
        if wildcards.coverage == 32.51:
            return ""
        else:
            return " -s " + seed + "." + str(float(wildcards.coverage) / 32.51)[2:]

def getScalingFactorOld(wildcards):
    seed = "862026"
    if wildcards.technology == "ont":
        if wildcards.coverage == 65.58:
            return ""
        else:
            return " -s " + seed + "." + str(float(wildcards.coverage) / 65.58)[2:]
    elif wildcards.technology == "pacbio":
        if wildcards.coverage == 35.35:
            return ""
        else:
            return " -s " + seed + "." + str(float(wildcards.coverage) / 35.35)[2:]

def getK(wildcards):
    if wildcards.technology == "ont":
        return 11
    else:
        return 21


# --- variable parameters within the main simulation ---#
coverages = [10, 15, 30]
aln_params = ["1_-1_-1_-1", "1_-3_-1_-1", "1_-2_-3_-1"]
aln_params_global = ["1_-1_-1_-1", "1_-3_-10_-1"]
technologies = ["ont", "pacbio"]
ont_accuracies = [0.98, 0.99]
extraction_modes = ["partial", "complete"]
sim_methods = ["ITER_POA", "abPOA", "KMEANS", "TRUE_HAP", "fuzzy", "ITER_POA_PERMISSIVE", "Molephase"]

# --- consensus sequence methods to use ---#
ont_coverages = [15, 20, 32.51, 94.1]
pacbio_coverages = [15, 20, 32.51]
real_methods = ["ITER_POA", "abPOA", "KMEANS", "fuzzy", "ITER_POA_PERMISSIVE", "Molephase"]
time_methods = ["ITER_POA", "abPOA", "KMEANS", "fuzzy"] 

# --- parameters for kmer size simulation
kmer_sizes = [5, 7, 11, 21, 31]
ruleorder: create_consensus_sequence_kmer_k > create_consensus_sequences

rule sim:
    input:
        expand(config["results_dir"] + "/simulation/ont/{accuracy}/{coverage}/{p}/evaluation_{method}_{mode}_{type}.tsv", method = sim_methods, accuracy = ont_accuracies, coverage = coverages, p = aln_params, mode = extraction_modes, type = "local"),
        expand(config["results_dir"] + "/simulation/pacbio/0.99/{coverage}/{p}/evaluation_{method}_{mode}_{type}.tsv", coverage = coverages, p = aln_params, method = sim_methods, mode = extraction_modes, type = "local"),
        expand(config["results_dir"] + "/simulation/{technology}/{accuracy}/{coverage}/{p}/evaluation_{method}_{mode}_{type}.tsv", method = sim_methods, accuracy = ont_accuracies, technology = technologies, coverage = coverages, p = aln_params_global, mode = "complete", type = "global"),
        expand(config["results_dir"] + "/simulation/{technology}/{accuracy}/{coverage}/1_-1_-1_-1/evaluation_{method}_{mode}_{type}.tsv", method = sim_methods, accuracy = 0.99, technology = technologies, coverage = coverages, mode = "partial", type = "semiglobal")

rule sim_giab:
    input:
        expand(config["results_dir"] + "/giab_simulation/{technology}/{coverage}/evaluation_{method}_{mode}_{type}.tsv", technology = technologies, method = sim_methods, coverage = [10, 15], mode = extraction_modes, type = "local"),
        expand(config["results_dir"] + "/giab_simulation/{technology}/{coverage}/evaluation_{method}_{mode}_{type}.tsv", technology = technologies, method = sim_methods, coverage = [10, 15], mode = "complete", type = "global"),
        expand(config["results_dir"] + "/giab_simulation/{technology}/{coverage}/evaluation_{method}_{mode}_{type}.tsv", technology = technologies, method = sim_methods, coverage = [10, 15], mode = "partial", type = "semiglobal")

rule time_benchmark:
    input:
        expand("data/time_benchmark/{mode}/{type}/consensus_stats_{method}_{technology}_{coverage}_{margin}.tsv", coverage = pacbio_coverages, method = time_methods, technology = technologies, mode = "partial", type = "local", margin = [500, 1000])

rule benchmark_k:
    input:
        expand(config["results_dir"] + "/kmeans_simulation/{technology}/0.99/{coverage}/{p}/evaluation_KMEANS_{k}_{mode}_{type}.tsv", technology = technologies, coverage = 15, p = "1_-1_-1_-1", k = kmer_sizes, mode = "partial", type = "local"),
        expand(config["results_dir"] + "/kmeans_simulation/{technology}/0.99/{coverage}/{p}/evaluation_KMEANS_{k}_{mode}_{type}.tsv", technology = technologies, coverage = 15, p = "1_-1_-1_-1", k = kmer_sizes, mode = "complete", type = "global"),
        expand(config["results_dir"] + "/kmeans_simulation/ont/0.98/{coverage}/{p}/evaluation_KMEANS_{k}_{mode}_{type}.tsv", coverage = 15, p = "1_-1_-1_-1", k = kmer_sizes, mode = "partial", type = "local"),
        expand(config["results_dir"] + "/kmeans_simulation/ont/0.98/{coverage}/{p}/evaluation_KMEANS_{k}_{mode}_{type}.tsv", coverage = 15, p = "1_-1_-1_-1", k = kmer_sizes, mode = "complete", type = "global")

rule real:
    input:
        expand(config["results_dir"] + "/real/pacbio/{coverage}/evaluation_{method}_partial_{mode}_{margin}.tsv", coverage = pacbio_coverages, method = real_methods, mode = ["local", "semiglobal"], margin = [500, 1000]),
        expand(config["results_dir"] + "/real/pacbio/{coverage}/evaluation_{method}_complete_{mode}_{margin}.tsv", coverage = pacbio_coverages, method = real_methods, mode = ["local", "global"], margin = [500, 1000]),
        expand(config["results_dir"] + "/real/ont/{coverage}/evaluation_{method}_partial_{mode}_{margin}.tsv", coverage = ont_coverages, method = real_methods, mode = ["local", "semiglobal"], margin = [500, 1000]),
        expand(config["results_dir"] + "/real/ont/{coverage}/evaluation_{method}_complete_{mode}_{margin}.tsv", coverage = ont_coverages, method = real_methods, mode = ["local", "global"], margin = [500, 1000])

rule real_old:
    input:
        expand(config["results_dir"] + "/real/pacbio_old/{coverage}/evaluation_{method}_partial_{mode}_{margin}.tsv", coverage = pacbio_coverages, method = real_methods, mode = ["local", "semiglobal"], margin = [500, 1000]),
        expand(config["results_dir"] + "/real/pacbio_old/{coverage}/evaluation_{method}_complete_{mode}_{margin}.tsv", coverage = pacbio_coverages, method = real_methods, mode = ["local", "global"], margin = [500, 1000]),
        expand(config["results_dir"] + "/real/ont_old/{coverage}/evaluation_{method}_partial_{mode}_{margin}.tsv", coverage = pacbio_coverages, method = real_methods, mode = ["local", "semiglobal"], margin = [500, 1000]),
        expand(config["results_dir"] + "/real/ont_old/{coverage}/evaluation_{method}_complete_{mode}_{margin}.tsv", coverage = pacbio_coverages, method = real_methods, mode = ["local", "global"], margin = [500, 1000])

rule plotRealCoverage:
    input:
        config["realdata_dir"] + "/t2t_remap/32.51/coverages.svg"

rule uncertainty:
    input:
        expand(config["results_dir"] + "/uncertainty_simulation/{technology}/{margin}/{jitter}/evaluation_{method}_{mode}_{type}.tsv", technology = technologies, margin = [100, 500, 1000], jitter = [0, 50, 200], method = sim_methods, mode = "partial", type = "local"), # margin 100
        expand(config["results_dir"] + "/uncertainty_simulation/{technology}/{margin}/{jitter}/evaluation_{method}_{mode}_{type}.tsv", technology = technologies, margin = [100, 500, 1000], jitter = [0, 50, 200], method = sim_methods, mode = "complete", type = "global") # margin 100


rule download_t2t:
    output:
        ref = config["real_ref"],
        fai = config["real_ref"] + ".fai"
    conda: "envs/minimap_env.yaml"
    log: "logs/download_t2t.log"
    shell:
        '(wget -O {output.ref}.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz && gunzip {output.ref}.gz && samtools faidx {output.ref}) &> {log}'

rule extract_chr1:
    input:
        ref = config["real_ref"]
    output:
        ref = config["ref_file"],
        fai = config["ref_file"] + ".fai"
    conda: "envs/minimap_env.yaml"
    log: "logs/extract_chr1.log"
    shell:
        '(samtools faidx {input.ref} chr1 > {output.ref} && samtools faidx {output.ref}) &> {log}'

include: "snakemake_modules/simulation.smk"
include: "snakemake_modules/real.smk"
#include: "snakemake_modules/real_old.smk"
