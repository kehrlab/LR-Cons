from intervaltree import IntervalTree
import sys


def add_standard_bed(repeat_tree, bedfile, reptype = "X"):
    with open(bedfile, 'r') as f_in:
        for line in f_in:
            linedata = line.strip().split('\t')
            chrom = linedata[0]
            start = int(linedata[1])
            end = int(linedata[2])
            if chrom not in repeat_tree:
                repeat_tree[chrom] = IntervalTree()
            repeat_tree[chrom][start:end] = reptype


if __name__ == "__main__":
    if len(sys.argv) != 9:
        print("Usage: python classify_regions.py <VARIANT_REGIONS.tsv> <RM.bed> <SD.full.bed> <CENSAT.bed> <TELOMERE.bed> <TRF.bed> <rDNA.bed> <OUTPUT.tsv>")
        sys.exit(-1)

    varfile = sys.argv[1]
    rm_regions = sys.argv[2]
    sd_regions = sys.argv[3]
    censat_regions = sys.argv[4]
    telomere_regions = sys.argv[5]
    trf_regions = sys.argv[6]
    rdna_regions = sys.argv[7]
    outfile = sys.argv[8]

    repeat_tree = {}

    # add RepeatMasker repeats to IntervalTrees
    # - SR
    # - ALU
    # - Any
    print("Processing RepeatMasker...")
    with open(rm_regions, 'r') as f_in:
        for line in f_in:
            linedata = line.strip().split('\t')
            chrom = linedata[0]
            start = int(linedata[1])
            end = int(linedata[2])
            type = linedata[6]
            subtype = linedata[7]

            if chrom not in repeat_tree:
                repeat_tree[chrom] = IntervalTree()

            cl = "RM"
            if type == "Simple_repeat":
                cl = "SR"
            elif type == "SINE" and subtype == "Alu":
                cl = "Alu"

            repeat_tree[chrom][start:end] = cl

    # add SDs
    print("Processing Segmental Duplications")
    with open(sd_regions, 'r') as f_in:
        for line in f_in:
            if line[0] == '#':
                continue
            linedata = line.strip().split('\t')
            chrom1 = linedata[0]
            start1 = int(linedata[1])
            end1 = int(linedata[2])
            chrom2 = linedata[9]
            start2 = int(linedata[10])
            end2 = int(linedata[11])

            if chrom1 not in repeat_tree:
                repeat_tree[chrom1] = IntervalTree()
            repeat_tree[chrom1][start1:end1] = "SD"

            if chrom2 not in repeat_tree:
                repeat_tree[chrom2] = IntervalTree()
            repeat_tree[chrom2][start2:end2] = "SD"

    # add censat regions
    print("Processing CenSat...")
    add_standard_bed(repeat_tree, censat_regions, 'C')
    print("Processing Telomeres...")
    add_standard_bed(repeat_tree, telomere_regions, 'T')
    print("Processing TRF...")
    add_standard_bed(repeat_tree, trf_regions, 'TR')
    print("Processing rDNA...")
    add_standard_bed(repeat_tree, rdna_regions, 'RD')


    print("Classify Variant Regions...")
    # classify variant regions based on Tree
    with open(outfile, 'w') as f_out:
        f_out.write("VARIANT\tNone\tRM\tSR\tAlu\tSD\tCenSat\tTelomere\tTRF\trDNA\n")
        with open(varfile, 'r') as f_in:
            for line in f_in:
                linedata = line.strip().split('\t')
                vname = linedata[0]
                regionstring = linedata[4][:-1]
                chrom = regionstring.split(':')[0]
                start = int(regionstring.split(':')[1].split('-')[0])
                end = int(regionstring.split(':')[1].split('-')[1])

                if chrom not in repeat_tree:
                    f_out.write(vname + "\t1\t0\t0\t0\t0\t0\t0\t0\t0\n")
                    continue
                tree_result = repeat_tree[chrom].overlap(start, end)
                if len(tree_result) == 0:
                    f_out.write(vname + "\t1\t0\t0\t0\t0\t0\t0\t0\t0\n")
                    continue

                rm_flag = 0
                sd_flag = 0
                sr_flag = 0
                alu_flag = 0
                censat_flag = 0
                tel_flag = 0
                rDNA_flag = 0
                trf_flag = 0

                for interval in tree_result:
                    overlap = (min(end, interval[1]) - max(start, interval[0])) / (end - start)
                    if interval[2] == 'RM':
                        rm_flag += overlap
                    elif interval[2] == 'SD':
                        sd_flag += overlap
                    elif interval[2] == 'SR':
                        sr_flag += overlap
                    elif interval[2] == 'Alu':
                        alu_flag += overlap
                        rm_flag += overlap
                    elif interval[2] == 'C':
                        censat_flag = 1
                    elif interval[2] == 'T':
                        tel_flag = 1
                    elif interval[2] == 'TR':
                        trf_flag += overlap
                    elif interval[2] == 'RD':
                        rDNA_flag = 1

                # if there is at least one Alu, check whether there is one overlapping either side of the repeat
                if alu_flag > 0:
                    alus = [(interval[0], interval[1]) for interval in tree_result if interval[2] == 'Alu']
                    left_overlap = [interval for interval in alus if interval[0] <= (start + 500) and interval[1] >= (start + 500)]
                    right_overlap = [interval for interval in alus if interval[0] <= (end - 500) and interval[1] >= (end - 500)]
                    if len(left_overlap) > 0 and len(right_overlap) > 0:
                        alu_flag = 1
                    else:
                        alu_flag = 0
                all_zero = rm_flag == 0 and sd_flag == 0 and sr_flag == 0 and alu_flag == 0

                # write classification
                f_out.write(
                    vname + '\t' +
                    str(int(all_zero)) + '\t' +
                    str(rm_flag) + '\t' +
                    str(sr_flag) + '\t' +
                    str(alu_flag) + '\t' +
                    str(sd_flag) + '\t' +
                    str(censat_flag) + '\t' +
                    str(tel_flag) + '\t' +
                    str(trf_flag) + '\t' +
                    str(rDNA_flag) + '\n'
                )
