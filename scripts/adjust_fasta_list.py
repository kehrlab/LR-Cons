#!/bin/bash

import sys

if len(sys.argv) != 3:
    exit(-1)

fa = sys.argv[1]
fout = sys.argv[2]
with open(fout, 'w') as f_out:
    with open(fa, 'r') as f_in:
        for line in f_in:
            data = line.strip().split('\t')
            region = data[2].split('_')
            if len(region) == 2 and region[0][:4] == 'chr1':
                region = 'chr1_' + region[0][4:] + '_' + region[1]
            f_out.write(data[0] + '\t' + data[1] + '\t' + region + '\n')
