configfile: "config.yaml"

# helper functions for parameter extraction
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

#--- definition of regons for Molephase ---#
def define_phasing_regions(size = 100000, overlap = 50000, chromosome = "chr1"):
    regions = []
    start = 0
    while start < 248010000:
        end = start + size
        regions.append(chromosome + "_" + str(start) + "_" + str(end))
        start += size - overlap
    return regions

#--- variable parameters within the benchmark ---#
coverages = [10, 15, 30]
aln_params = ["1_-1_-1_-1", "1_-3_-1_-1", "1_-2_-3_-1"]
aln_params_global = ["1_-1_-1_-1", "1_-3_-10_-1"]
technologies = ["ont", "pacbio"]
ont_accuracies = [0.99]
#--- consensus sequence methods to use ---#
methods = ["ITER_POA", "abPOA", "KMEANS", "Molephase", "TRUE_HAP", "fuzzy", "ITER_POA_PERMISSIVE"]

rule all:
    input:
        expand(config["results_dir"] + "/ont/{accuracy}/{coverage}/{p}/evaluation_{method}.tsv", method = methods, accuracy = ont_accuracies, coverage = coverages, p = aln_params),
        expand(config["results_dir"] + "/pacbio/0.99/{coverage}/{p}/evaluation_{method}.tsv", coverage = coverages, p = aln_params, method = methods),
        expand(config["results_dir"] + "/{technology}/{accuracy}/{coverage}/{p}/evaluation_{method}_global.tsv", method = methods, accuracy = ont_accuracies, technology = technologies, coverage = coverages, p = aln_params_global),
        expand(config["results_dir"] + "/{technology}/{accuracy}/{coverage}/1_-1_-1_-1/evaluation_{method}_semiglobal.tsv", method = methods, accuracy = ont_accuracies, technology = technologies, coverage = coverages),
        expand(config["results_dir"] + "/{technology}/{accuracy}/{coverage}/1_-1_-1_-1/evaluation_{method}_spanning.tsv", method = methods, accuracy = ont_accuracies, technology = technologies, coverage = coverages)


rule time_benchmark:
    input:
        expand(config["data_dir"] + "/time_data/{coverage}/consensus_stats_{method}_{technology}_{accuracy}_{p}.tsv", coverage = coverages, method = methods, technology = technologies, accuracy = ont_accuracies, p = aln_params)

include: "snakemake_fragments/variant_sim.smk"
