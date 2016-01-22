#!/usr/bin/env python

import csv
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--definitions_file', help="Illumina Definitions File. Variants Only")
    parser.add_argument('-o', '--output', help="Root name for output files. Separate files for SNVs/Substitutions"
                                               " and Indels created")
    args = parser.parse_args()

    simple_variants = dict()
    indels = dict()

    with open(args.definitions_file, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:
            if len(row[3]) == len(row[4]):
                simple_variants["{chr}\t{start}\t{stop}".format(chr=row[1],
                                                                start=(int(row[2]) - 1),
                                                                stop=(int(row[2]) + len(row[3])))] = row[0]
            else:
                if len(row[3]) > len(row[4]):
                    size = len(row[3])
                else:
                    size = len(row[4])

                indels["{chr}\t{start}\t{stop}".format(chr=row[1], start=row[2], stop=(int(row[2]) + size))] = row[0]

    with open("{}.indel.bed".format(args.output), 'w') as indels_out:
        for region in indels:
            indels_out.write("{region}\t{gene}\n".format(region=region, gene=indels[region]))

    with open("{}.snv.bed".format(args.output), 'w') as snvs_out:
        for region in simple_variants:
            snvs_out.write("{region}\t{gene}\n".format(region=region, gene=simple_variants[region]))
