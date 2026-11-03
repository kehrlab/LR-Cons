import multiprocessing
import os
import subprocess
import sys
import pysam


def get_chromosome_info(ref_path):
    chromosome_lens = {}
    with open(ref_path + ".fai", "r") as f_in:
        for line in f_in:
            values = line.strip().split("\t")
            chromosome_lens[values[0]] = int(values[1])
    return chromosome_lens


# --- definition of regons for Molephase ---#
def define_phasing_regions(cname, clen, size=100000, overlap=50000):
    regions = []
    start = 0
    while start < clen:
        end = start + size
        regions.append({"cname": cname, "start": start, "end": end})
        start += size - overlap
    return regions


def phase_chromosome(chrom, clen, molephase_bin, bam_path, ref_path, out_path, k):
    regions = define_phasing_regions(chrom, clen, 100000, 50000)
    for region in regions:
        region_string = (
            region["cname"] + ":" + str(region["start"]) + "-" + str(region["end"])
        )
        prefix = (
            out_path
            + "/"
            + region["cname"]
            + "_"
            + str(region["start"])
            + "_"
            + str(region["end"])
        )

        all_path = prefix + "_all.fq"

        cmd_string = (
            "samtools view -h "
            + bam_path
            + " "
            + region_string
            + " | samtools fastq - > "
            + all_path
            + "; ./"
            + molephase_bin
            + " -1 "
            + all_path
            + " -o "
            + prefix
            + " -r "
            + ref_path
            + " -c "
            + region["cname"]
            + " -s "
            + str(region["start"])
            + " -e "
            + str(region["end"])
            + " -k "
            + str(k)
        )
        os.system(cmd_string)
        try:
            pysam.faidx(prefix + ".hap1.fq")
        except:
            print("Could not index " + prefix + ".hap1.fq")
        try:
            pysam.faidx(prefix + ".hap2.fq")
        except:
            print("Could not index " + prefix + ".hap2.fq")
        try:
            pysam.faidx(prefix + ".unphased.fq")
        except:
            print("Could not index " + prefix + ".unphased.fq")

    return True


if __name__ == "__main__":
    molephase_bin = sys.argv[1]
    bam_path = sys.argv[2]
    ref_path = sys.argv[3]
    out_path = sys.argv[4]
    k = int(sys.argv[5])
    threads = int(sys.argv[6])

    cinfo = get_chromosome_info(ref_path)

    os.mkdir(out_path)

    items = [
        (chrom, cinfo[chrom], molephase_bin, bam_path, ref_path, out_path, k)
        for chrom in cinfo
    ]
    with multiprocessing.Pool(processes=threads) as pool:
        pool.starmap(phase_chromosome, items)
