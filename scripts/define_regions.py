import json
import sys

if len(sys.argv) != 4:
    print("Invalid call: ./define_regions.py <VARIANT_INFO> <MARGIN> <OUT_FILE> ")
    sys.exit(-1)

info_json = sys.argv[1]
margin = int(sys.argv[2])
out_file = sys.argv[3]

with open(info_json, "r") as f_in:
    variants = json.load(f_in)

breakpoints = {}
sizes = {}
types = {}
zygosities = {}
chromosomes = {}

for id in variants:
    variant = variants[id]
    if ".1" in id or ".2" in id:
        insert_zyg = ".1" in id
        id = id[:-2]
        if insert_zyg:
            zygosities[id] = variant["ZYG"]
    else:
        zygosities[id] = variant["ZYG"]

    bp1 = variant["START"]
    bp2 = variant["END"]

    if id not in breakpoints:
        breakpoints[id] = [bp1, bp2]
    else:
        breakpoints[id] += [bp1, bp2]

    # zygosities[id] = variant["ZYG"]
    types[id] = variant["TYPE"]
    chromosomes[id] = variant["CHR"]

# convert breakpoint list to regions and merge overlapping regions
regions = {}
for id in breakpoints:
    regions[id] = []
    breakpoints[id].sort()
    region_start = breakpoints[id][0] - margin
    region_end = breakpoints[id][0] + margin
    for i in range(len(breakpoints[id])):
        region_end = breakpoints[id][i] + margin
        if i == len(breakpoints[id]) - 1:
            break
        if region_end < breakpoints[id][i + 1] - margin:
            regions[id].append("chr1:" + str(region_start) + "-" + str(region_end))
            region_start = breakpoints[id][i + 1] - margin
    regions[id].append("chr1:" + str(region_start) + "-" + str(region_end))

# write variant info and regions
with open(out_file, "w") as f_out:
    for id in breakpoints:
        f_out.write(
            id
            + "\t"
            + chromosomes[id]
            + "\t"
            + types[id]
            + "\t"
            + zygosities[id]
            + "\t"
        )
        # for bp in breakpoints[id]:
        #     f_out.write("\t" + str(bp))
        for region in regions[id]:
            f_out.write(region + ",")
        f_out.write("\n")
