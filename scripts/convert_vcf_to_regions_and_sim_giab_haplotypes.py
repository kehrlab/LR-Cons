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

    def getREF(self):
        return self.ref

    def getSeq(self):
        return self.alt


def insert_snvs(ref_seqs, sm_infile):
    haplotype_seqs = {}

    current_chr = "void"
    hap1_pos = 0
    hap2_pos = 0
    hap1_snippets = []
    hap2_snippets = []
    with open(sm_infile, 'r') as f_in:
        for line in f_in:
            # skip the header
            if line[0] == "#":
                continue
            record = VcfRecord(line.split())
            chrom = record.chrom
            if chrom != current_chr:
                # need to complete the previous chromosome sequence
                if current_chr != "void":
                    hap1_snippets.append(ref_seqs[current_chr][hap1_pos:])
                    hap2_snippets.append(ref_seqs[current_chr][hap2_pos:])
                    haplotype_seqs[current_chr][0] = ''.join(hap1_snippets)
                    haplotype_seqs[current_chr][1] = ''.join(hap2_snippets)
                #
                hap1_snippets = []
                hap2_snippets = []
                hap1_pos = 0
                hap2_pos = 0
                haplotype_seqs[chrom] = ["", ""]
                current_chr = chrom
            p = record.pos

            gt = record.getGT()
            if '.' in gt or '|' not in gt:
                continue

            refseq = record.getREF()
            seq = record.getSeq()
            if len(refseq) != 1:
                continue
            valid_alt = True
            for s in seq:
                if len(s) != 1:
                    valid_alt = False
            if not valid_alt:
                continue

            # add the sequence for each haplotype
            if gt[0] != "0":
                alt_idx = int(gt[0]) - 1
                hap1_snippets.append(ref_seqs[chrom][hap1_pos:(p-1)] + seq[alt_idx])
                # haplotype_seqs[chrom][0].join(ref_seqs[chrom][hap1_pos:(p-1)], seq[alt_idx])
                hap1_pos = p
            if gt[2] != "0":
                alt_idx = int(gt[2]) - 1
                hap2_snippets.append(ref_seqs[chrom][hap2_pos:(p-1)] + seq[alt_idx])
                # haplotype_seqs[chrom][1] = haplotype_seqs[chrom][1] + ref_seqs[chrom][hap2_pos:(p-1)] + seq[alt_idx]
                hap2_pos = p
        hap1_snippets.append(ref_seqs[current_chr][hap1_pos:])
        hap2_snippets.append(ref_seqs[current_chr][hap2_pos:])
        haplotype_seqs[chrom][0] = ''.join(hap1_snippets)
        haplotype_seqs[chrom][1] = ''.join(hap2_snippets)
    return haplotype_seqs


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
    if type == "DEL":
        sv_seq = ref_seqs[chrom][max(0, insert_pos - margin) : (insert_pos)]
        sv_seq += ref_seqs[chrom][(insert_pos + l) : (insert_pos + l + margin)]
        allele_seqs.append(sv_seq)
    elif type == "INS":
        sv_seq = ref_seqs[chrom][max(0, insert_pos - margin) : (insert_pos - 1)]
        temp_seq = ref_seqs[chrom][(insert_pos) : (insert_pos + margin)]
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

    svs = dict()
    for k in ref_seqs:
        svs[k] = []

    with open(inFile, "r") as fIn:
        for line in fIn:
            # skip the header
            if line[0] == "#":
                continue

            record = VcfRecord(line.split())
            chrom = record.chrom

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
            # alleles = getAlleleSeqs(record, ref_seqs, 4 * margin)

            sv_details = {
                "pos" : p,
                "type": record.info["SVTYPE"],
                "size" : svlen,
                "seq": record.getSeq().split(';'),
                "hap1": gt[0],
                "hap2": gt[2],
                "idx": idx
            }
            svs[chrom].append(sv_details)

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
            # sequences[id] = alleles
    return svs #sequences, svs


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print(
            "Usage: ./convert_vcf_to_regions.py <VCF_IN> <SNV_VCF_IN> <REFERENCE> <MARGIN> <ALLELE_TRUTH_DIR> <HAPLOTYPE_FILE>"
        )
        sys.exit()

    infile = sys.argv[1]
    sm_infile = sys.argv[2]
    reference = sys.argv[3]
    margin = int(sys.argv[4])
    # seq_file = sys.argv[4]
    truth_dir = sys.argv[5]
    out_file = sys.argv[6]
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

    svs = convertVCFfile(infile, ref_seqs, margin)


    # create full haplotype sequences...
    haplotype_seqs = insert_snvs(ref_seqs, sm_infile)


    sv_positions = {}

    for chrom in haplotype_seqs:
        hap1_snippets = []
        hap2_snippets = []
        hap1_pos = 0
        hap2_pos = 0

        hap1_offset = 0
        hap2_offset = 0


        for sv_details in svs[chrom]:
            id = chrom + "_" + sv_details["type"] + "_" + str(sv_details['pos']) + "_" + str(sv_details['size']) + "_" + str(sv_details['idx'])
            sv_positions[id] = {
                "h1_pos" : sv_details["pos"] + hap1_offset,
                "h1_end" : sv_details["pos"] + hap1_offset + 1,
                "h2_pos" : sv_details["pos"] + hap2_offset,
                "h2_end" : sv_details["pos"] + hap2_offset + 1,
                "size": sv_details["size"],
                "type": sv_details["type"],
                "chrom": chrom
            }
            if sv_details["type"] == "DEL":
                if sv_details['hap1'] != "0":
                    hap1_snippets.append(haplotype_seqs[chrom][0][hap1_pos:sv_details['pos']])
                    hap1_pos = sv_details['pos'] + abs(sv_details['size'])
                    hap1_offset -= abs(sv_details['size'])
                if sv_details['hap2'] != "0":
                    hap2_snippets.append(haplotype_seqs[chrom][1][hap2_pos:sv_details['pos']])
                    hap2_pos = sv_details['pos'] + abs(sv_details['size'])
                    hap2_offset -= abs(sv_details['size'])
            else:
                if sv_details['hap1'] != "0":
                    sv_positions[id]["h1_end"] += abs(sv_details["size"])
                    hap1_snippets.append(haplotype_seqs[chrom][0][hap1_pos:(sv_details['pos']-1)])
                    hap1_snippets.append(sv_details['seq'][int(sv_details['hap1']) - 1])
                    hap1_pos = sv_details['pos']
                    hap1_offset += abs(sv_details['size'])
                if sv_details['hap2'] != "0":
                    sv_positions[id]["h2_end"] += abs(sv_details["size"])
                    hap2_snippets.append(haplotype_seqs[chrom][1][hap2_pos:(sv_details['pos']-1)])
                    hap2_snippets.append(sv_details['seq'][int(sv_details['hap2']) - 1])
                    hap2_pos = sv_details['pos']
                    hap2_offset += abs(sv_details['size'])
        hap1_snippets.append(haplotype_seqs[chrom][0][hap1_pos:])
        hap2_snippets.append(haplotype_seqs[chrom][1][hap2_pos:])
        haplotype_seqs[chrom][0] = ''.join(hap1_snippets)
        haplotype_seqs[chrom][1] = ''.join(hap2_snippets)


    for id in sv_positions:
        chrom = sv_positions[id]["chrom"]
        hap_seqs = [
            haplotype_seqs[chrom][0][max(0, sv_positions[id]["h1_pos"] - 2000) : (sv_positions[id]["h1_end"] + 2000)],
            haplotype_seqs[chrom][1][max(0, sv_positions[id]["h2_pos"] - 2000) : (sv_positions[id]["h2_end"] + 2000)]
        ]

        truth_file = truth_dir + "/" + id + "_truth.fa"
        with open(truth_file, "w") as f_out:
            for i in range(len(hap_seqs)):
                s = hap_seqs[i]
                f_out.write("> HAP_" + str(i) + "\n")
                p = 0
                l = 70
                while p < len(s):
                    f_out.write(s[p : p + l] + "\n")
                    p += l
        pysam.faidx(truth_file)

    if out_file != "/dev/null":
        # write haplotype sequences
        with open(out_file, 'w') as f_out:
            for chrom in haplotype_seqs:
                for i in range(2):
                    id = chrom + ".hap_" + str(i)
                    f_out.write("> " + id + "\n")
                    p = 0
                    l = 70
                    while p < len(haplotype_seqs[chrom][i]):
                        f_out.write(haplotype_seqs[chrom][i][p : p+l] + "\n")
                        p += l
        pysam.faidx(out_file)
