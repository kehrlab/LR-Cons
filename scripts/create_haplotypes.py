import json
import os
import random
import sys

import pysam
from intervaltree import IntervalTree
from tqdm import tqdm

comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def reverse_complement(seq):
    return "".join([comp[s] for s in seq[::-1]])


def insert_deletion(haplotype_seqs, chromosome, size, i, start, end):
    # determine zygosity
    hap = random.choice(
        [1, 2, 3]
    )  # 1: first haplotype, 2: second, 3: both -> 2:1 ratio of heterozygous to homozygous

    id = "DEL_" + str(i) + "_" + str(size)
    variants = {
        id: {
            "CHR": chromosome,
            "START": start,
            "END": end,
            "SIZE": size,
            "ZYG": ["", "HAP1", "HAP2", "HOM"][hap],
            "TYPE": "DEL",
        }
    }

    for i in range(len(haplotype_seqs)):
        if (i == 0 and hap != 2) or (i == 1 and hap != 1):
            haplotype_seqs[i] = (
                haplotype_seqs[i][: start + 1] + haplotype_seqs[i][end + 1 :]
            )
    return variants


def insert_inversion(haplotype_seqs, chromosome, size, i, start, end):
    # determine zygosity
    hap = random.choice(
        [1, 2, 3]
    )  # 1: first haplotype, 2: second, 3: both -> 2:1 ratio of heterozygous to homozygous

    id = "INV_" + str(i) + "_" + str(size)
    variants = {
        id: {
            "CHR": chromosome,
            "START": start,
            "END": end,
            "SIZE": size,
            "ZYG": ["", "HAP1", "HAP2", "HOM"][hap],
            "TYPE": "INV",
        }
    }

    for i in range(len(haplotype_seqs)):
        if (i == 0 and hap != 2) or (i == 1 and hap != 1):
            haplotype_seqs[i] = (
                haplotype_seqs[i][: start + 1]
                + reverse_complement(haplotype_seqs[i][start + 1 : end + 1])
                + haplotype_seqs[i][end + 1 :]
            )
    return variants


def insert_duplication(haplotype_seqs, chromosome, size, i, start, end):
    # determine zygosity
    hap = random.choice(
        [1, 2, 3]
    )  # 1: first haplotype, 2: second, 3: both -> 2:1 ratio of heterozygous to homozygous

    id = "DUP_" + str(i) + "_" + str(size)
    variants = {
        id: {
            "CHR": chromosome,
            "START": start,
            "END": end,
            "SIZE": size,
            "ZYG": ["", "HAP1", "HAP2", "HOM"][hap],
            "TYPE": "DUP",
        }
    }

    for i in range(len(haplotype_seqs)):
        if (i == 0 and hap != 2) or (i == 1 and hap != 1):
            haplotype_seqs[i] = (
                haplotype_seqs[i][: end + 1] + haplotype_seqs[i][start + 1 :]
            )
    return variants


def insert_comphet_del(haplotype_seqs, chromosome, size, i, start, end):
    id1 = "COMPHET_DEL_" + str(i) + "_" + str(size) + ".1"
    id2 = "COMPHET_DEL_" + str(i) + "_" + str(size) + ".2"
    variants = {
        id1: {
            "CHR": chromosome,
            "START": start,
            "END": end,
            "SIZE": size,
            "ZYG": "",
            "TYPE": "COMPHET_DEL",
        },
        id2: {
            "CHR": chromosome,
            "START": start - int(0.1 * size),
            "END": end + int(0.2 * size),
            "SIZE": size,
            "ZYG": "",
            "TYPE": "COMPHET_DEL",
        },
    }

    # create one deletion according to specification and one that is slightly shifted and 30% longer
    indices = random.sample([0, 1], k=2)
    haplotype_seqs[indices[0]] = (
        haplotype_seqs[indices[0]][: start + 1] + haplotype_seqs[indices[0]][end + 1 :]
    )
    variants[id1]["ZYG"] = ["HAP1", "HAP2"][indices[0]]

    haplotype_seqs[indices[1]] = (
        haplotype_seqs[indices[1]][: start + 1 - int(0.1 * size)]
        + haplotype_seqs[indices[1]][
            start + size + int(0.3 * size) + 1 - int(0.1 * size) :
        ]
    )
    variants[id2]["ZYG"] = ["HAP1", "HAP2"][indices[1]]

    return variants


def insert_shifted_comphet(haplotype_seqs, chromosome, size, i, start, end):
    id1 = "SHIFTED_COMPHET_" + str(i) + "_" + str(size) + ".1"
    id2 = "SHIFTED_COMPHET_" + str(i) + "_" + str(size) + ".2"
    variants = {
        id1: {
            "CHR": chromosome,
            "START": start,
            "END": end,
            "SIZE": size,
            "ZYG": "",
            "TYPE": "SHIFTED_COMPHET",
        },
        id2: {
            "CHR": chromosome,
            "START": end + 150,
            "END": end + 150 + int(size * 0.6),
            "SIZE": int(size * 0.6),
            "ZYG": "",
            "TYPE": "SHIFTED_COMPHET",
        },
    }

    indices = random.sample([0, 1], k=2)

    haplotype_seqs[indices[0]] = (
        haplotype_seqs[indices[0]][: start + 1] + haplotype_seqs[indices[0]][end + 1 :]
    )
    variants[id1]["ZYG"] = ["HAP1", "HAP2"][indices[0]]

    haplotype_seqs[indices[1]] = (
        haplotype_seqs[indices[1]][: end + 1 + 150]
        + reverse_complement(
            haplotype_seqs[indices[1]][end + 1 + 150 : end + 1 + 150 + int(size * 0.6)]
        )
        + haplotype_seqs[indices[1]][end + 1 + 150 + int(size * 0.6) :]
    )
    variants[id2]["ZYG"] = ["HAP1", "HAP2"][indices[1]]

    return variants


def get_allele_seqs(haplotype_seqs, variant_info, margin=1000):
    start = 250000000
    end = 0

    allele_seqs = ["", ""]

    comphet = False
    for v in variant_info:
        start = min(start, variant_info[v]["START"])
        if (
            variant_info[v]["TYPE"] == "COMPHET_DEL"
            or variant_info[v]["TYPE"] == "SHIFTED_COMPHET"
        ):
            comphet = True

    if comphet:
        for v in variant_info:
            if variant_info[v]["TYPE"] == "SHIFTED_COMPHET":
                if ".1" in v:
                    end = start + int(0.6 * variant_info[v]["SIZE"]) + 150
                else:
                    end = start + int(variant_info[v]["SIZE"] * (1.6 / 0.6)) + 150
            else:
                if ".1" in v:
                    end = start + int(0.2 * variant_info[v]["SIZE"])
                else:
                    end = start + 1

            hap_idx = 0
            if variant_info[v]["ZYG"] == "HAP2":
                hap_idx = 1
            allele_seqs[hap_idx] = haplotype_seqs[hap_idx][
                (start - margin) : (end + margin)
            ]
    else:
        v = list(variant_info.keys())[0]
        ref_end = start + variant_info[v]["SIZE"]
        hap_end = start
        if variant_info[v]["TYPE"] == "INV":
            hap_end = ref_end
        elif variant_info[v]["TYPE"] == "DUP":
            hap_end = ref_end + variant_info[v]["SIZE"]
        start -= margin
        ref_end += margin
        hap_end += margin

        allele_seqs = [h[start:hap_end] for h in haplotype_seqs]
        if variant_info[v]["ZYG"] == "HAP1":
            allele_seqs[1] = haplotype_seqs[1][start:ref_end]
        elif variant_info[v]["ZYG"] == "HAP2":
            allele_seqs[0] = haplotype_seqs[0][start:ref_end]

    return allele_seqs


def create_snps(haplotype_seqs, e=0.0005):
    replacement_bases = {
        "A": ["C", "T", "G"],
        "C": ["A", "T", "G"],
        "T": ["A", "C", "G"],
        "G": ["A", "C", "T"],
    }
    for i in range(len(haplotype_seqs)):
        pbar = tqdm(range(len(haplotype_seqs[i])))
        pbar.set_description("Inserting SNPs in Haplotype " + str(i))
        temp_string = ""
        last_pos = 0
        n_snp = 0
        for j in pbar:
            if random.uniform(0, 1) <= e:
                temp_string += haplotype_seqs[i][last_pos:j] + random.choice(
                    replacement_bases[haplotype_seqs[i][j]]
                )
                last_pos = j + 1
                n_snp += 1
        temp_string += haplotype_seqs[i][last_pos:]
        haplotype_seqs[i] = temp_string
        print("\tInserted " + str(n_snp) + " SNPs")


def write_fasta(sequences, names, filename):
    with open(filename, "w") as f_out:
        for i in range(2):
            f_out.write(">" + names[i] + "\n")
            p = 0
            w = 70
            while p + 70 < len(sequences[i]):
                f_out.write(sequences[i][p : p + w] + "\n")
                p += w
            if p != len(sequences[i]):
                f_out.write(sequences[i][p:] + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print(
            "Invalid number of arguments: <REF> <EXCLUDE_REGIONS> <CHR> <INFO_OUT> <FASTA_OUT> <TRUTH_DIR>"
        )
        sys.exit(-1)

    random.seed(12341234)

    ref_file = sys.argv[1]
    exclude_bed = sys.argv[2]
    chromosome = sys.argv[3]
    info_file = sys.argv[4]
    fasta_file = sys.argv[5]
    truth_dir = sys.argv[6]

    if not os.path.exists(truth_dir):
        os.mkdir(truth_dir)
    else:
        if not os.path.isdir(truth_dir):
            print(truth_dir + " exists but does not appear to be a directoy!")
            sys.exit(-1)

    # read reference sequence
    ff = pysam.FastaFile(ref_file)
    if not ff.is_open():
        print("Could not open FASTA file ", ref_file)
        sys.exit(-1)

    ref_seqs = [ff.fetch(chromosome).upper()]
    create_snps(ref_seqs, 0.0003656)  # insert homozygous SNPs randomly
    print("Loaded ", len(ref_seqs), " reference sequences.")

    # read regions to exclude
    exclude_regions = []
    with open(exclude_bed, "r") as bed_in:
        for line in bed_in:
            data = line.strip().split()
            if data[0] != chromosome:
                continue
            positions = [int(p) for p in data[1:]]
            exclude_regions.append(positions)

    # fill IntervalTree with banned regions
    tree = IntervalTree()
    for r in exclude_regions:
        tree[r[0] : r[1]] = "X"

    # create variants
    variants = {}
    haplotype_seqs = [ref_seqs[0], ref_seqs[0]]
    start_pos = 248000000  # last unchanged positon; end: last changed base
    step_size = 15000

    variant_functions = [
        insert_deletion,
        insert_inversion,
        insert_duplication,
        insert_comphet_del,
        insert_shifted_comphet,
    ]

    n_cycles = 800  # 500
    sizes = [50, 200, 500, 1000]  # [50, 200, 500, 1000]
    # cycle through every variant type and size and insert variants
    for i in tqdm(range(n_cycles)):
        # print("Cycle ", i)
        for size in sizes:
            # print("\tSize ", size)
            for f in variant_functions:
                # print("\t\t", f.__name__)
                end_pos = start_pos + size
                # check [start_pos, start_pos + size] for overlap with exluded regions
                while len(tree.overlap(start_pos, end_pos)) > 0:
                    start_pos -= step_size
                    end_pos = start_pos + size

                if start_pos < 0:
                    print(
                        "Invalid position (< 0). Reduce step size or number of variants."
                    )
                    sys.exit(-1)

                variant_info = f(
                    haplotype_seqs, chromosome, size, i, start_pos, end_pos
                )
                allele_seqs = get_allele_seqs(haplotype_seqs, variant_info, 1000)

                type = ""
                for v in variant_info:
                    variants[v] = variant_info[v]
                    type = variants[v]["TYPE"]

                write_fasta(
                    allele_seqs,
                    ["HAP_1", "HAP_2"],
                    truth_dir
                    + "/"
                    + type
                    + "_"
                    + str(i)
                    + "_"
                    + str(size)
                    + "_truth.fa",
                )

                start_pos -= step_size

    # create heterozygous SNPs independently on the two alleles
    create_snps(haplotype_seqs, 0.0003172)

    # write results
    with open(info_file, "w") as f_out:
        json.dump(variants, f_out)

    hap_names = [chromosome + "_hap1", chromosome + "_hap2"]
    write_fasta(haplotype_seqs, hap_names, fasta_file)
