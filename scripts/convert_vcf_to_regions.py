import os
import sys

import pysam


class VcfRecord:
    def __init__(self, fields):
        self.chrom = fields[0]
        self.pos = int(fields[1])
        # self.end = self.pos
        self.id = fields[2]
        self.ref = fields[3]
        self.alt = fields[4]  ##only necessary for GRIDSS
        self.found = []

        qual = fields[5]
        if qual == ".":
            self.qual = None
        else:
            self.qual = round(float(qual))

        info = fields[7].split(";")
        self.info = {}
        for i in info:
            i = i.split("=")
            if len(i) == 1:
                i += [True]
            self.info[i[0]] = i[1]

        self.gt = fields[9].split(":")[0]

    def getChrom(self):
        return self.chrom

    def getPos(self):
        return self.pos

    def getLength(self):
        return abs(int(self.info["SVLEN"]))

    def isSV(self):
        return "SVTYPE" in self.info

    def getGT(self):
        return self.gt

    def getType(self):
        return self.info["SVTYPE"]

    def getSeq(self):
        return self.alt


def getAlleleSeqs(record, ref_seqs, margin):
    allele_seqs = []

    chrom = record.getChrom()
    insert_pos = record.getPos()
    alt = record.getSeq().split(";")
    type = record.getType()
    l = record.getLength()

    # reference sequence
    allele_seqs.append(
        ref_seqs[chrom][max(0, insert_pos - margin) : (insert_pos + margin)]
    )

    # allele_sequence
    sv_seq = ref_seqs[chrom][max(0, insert_pos - margin) : (insert_pos + 1)]
    if type == "DEL":
        sv_seq += ref_seqs[chrom][(insert_pos + l + 1) : (insert_pos + l + margin + 1)]
        allele_seqs.append(sv_seq)
    elif type == "INS":
        temp_seq = ref_seqs[chrom][(insert_pos + 1) : (insert_pos + 1 + margin)]
        for a in alt:
            allele_seqs.append(sv_seq + a + temp_seq)

    return allele_seqs


def convertVCFfile(inFile, ref_seqs, margin):
    """
    open input VCF
    """
    sequences = {}
    idx = 0

    too_small = 0
    too_large = 0
    no_gt = 0
    no_gt_trf = 0
    no_gt_rm = 0
    no_gt_easy = 0
    kept = 0
    kept_trf = 0
    kept_rm = 0
    kept_easy = 0

    with open(inFile, "r") as fIn:
        for line in fIn:
            # skip the header
            if line[0] == "#":
                continue

            record = VcfRecord(line.split())

            # keep only sufficiently large SVs
            if not record.isSV():
                continue

            svlen = record.getLength()
            if not svlen >= 50:
                too_small += 1
                continue

            if not svlen <= 5000:
                too_large += 1
                continue

            gt = record.getGT()
            if "." in gt or "|" not in gt:
                no_gt += 1
                if "TRF" in record.info:
                    no_gt_trf += 1
                if "RM_clsfam" in record.info:
                    no_gt_rm += 1
                if "TRF" not in record.info and "RM_clsfam" not in record.info:
                    no_gt_easy += 1
                continue

            chrom = record.chrom
            p = record.pos
            breakpoints = [p]

            if record.info["SVTYPE"] == "DEL":
                breakpoints.append(p + svlen)
            elif record.info["SVTYPE"] == "INS":
                breakpoints.append(p + 1)
            else:
                continue

            kept += 1
            if "TRF" in record.info:
                kept_trf += 1
            if "RM_clsfam" in record.info:
                kept_rm += 1
            if "TRF" not in record.info and "RM_clsfam" not in record.info:
                kept_easy += 1

            idx = idx + 1
            region = [max(breakpoints[0] - margin, 0), breakpoints[1] + margin]
            id = chrom + "_" + record.info["SVTYPE"] + "_" + str(p) + "_" + str(svlen) + "_" + str(idx)
            alleles = getAlleleSeqs(record, ref_seqs, 4 * margin)

            print(
                id
                + "\t"
                + chrom
                + "\t"
                + record.info["SVTYPE"]
                + "\t"
                + record.getGT()
                + "\t"
                + chrom
                + ":"
                + str(region[0])
                + "-"
                + str(region[1])
                + ","
            )
            sequences[id] = alleles
        # print(
        #     "Retained variants: ",
        #     kept,
        #     "; TRF: ",
        #     kept_trf,
        #     ", RM: ",
        #     kept_rm,
        #     ", other: ",
        #     kept_easy,
        # )
        # print("Too small: ", too_small)
        # print("Too large: ", too_large)
        # print(
        #     "Incomplete Genotype: ",
        #     no_gt,
        #     "; TRF: ",
        #     no_gt_trf,
        #     ", RM: ",
        #     no_gt_rm,
        #     ", other: ",
        #     no_gt_easy,
        # )
    return sequences


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(
            "Usage: ./convert_vcf_to_regions.py <VCF_IN> <REFERENCE> <MARGIN> <ALLELE_TRUTH_DIR>"
        )
        sys.exit()

    infile = sys.argv[1]
    reference = sys.argv[2]
    margin = int(sys.argv[3])
    # seq_file = sys.argv[4]
    truth_dir = sys.argv[4]
    if not os.path.exists(truth_dir):
        os.mkdir(truth_dir)

    if not os.path.exists(truth_dir):
        print("Could not create output directory: ", truth_dir)
        sys.exit(-1)

    ff = pysam.FastaFile(reference)
    if not ff.is_open():
        print("Could not open FASTA file ", reference)
        sys.exit(-1)

    ref_seqs = {}
    for chrom in ff.references:
        ref_seqs[chrom] = ff.fetch(chrom).upper()

    sequences = convertVCFfile(infile, ref_seqs, margin)

    # with open(seq_file, "w") as f_out:
    for id in sequences:
        if len(sequences[id]) > 2:
            print("Too many allele sequences. Not prepared to handle more than two.")
            sys.exit(-1)
        truth_file = truth_dir + "/" + id + "_truth.fa"
        outstring = id
        with open(truth_file, "w") as f_out:
            for i in range(len(sequences[id])):
                s = sequences[id][i]
                f_out.write("> HAP_" + str(i) + "\n")
                p = 0
                l = 70
                while p < len(s):
                    f_out.write(s[p : p + l] + "\n")
                    p += l
