#!/usr/bin/env python

import os
import argparse
import pybedtools

from collections import defaultdict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--depth_output', help="Output file name for depth information")
    parser.add_argument('-f', '--failed_output', help="Output file name for regions with failed samples")
    args = parser.parse_args()

    failed_regions = defaultdict(list)
    depth_regions = defaultdict(lambda: defaultdict(lambda: 'N/A'))
    samples_list = list()

    files = os.listdir(os.getcwd())
    for coverage_file in files:
        sections = coverage_file.partition(".")
        samples_list.append(sections[0])
        if sections[2] == 'diagnosetargets.vcf':
            coverage_data = pybedtools.BedTool(coverage_file)
            for region in coverage_data:
                info = region[7].split(';')
                data = region[9].split(":")
                filter_field = region[6]

                end = info[0].split('=')
                depth = data[-3]

                region_name = "{chr}:{start}-{end}".format(chr=region.chrom, start=region.start, end=end[1])
                depth_regions[region_name][sections[0]] = depth

                if filter_field != 'PASS':
                    failed_regions[region_name].append(sections[0])

    with open(args.failed_output, 'w') as failed_output:
        for region in failed_regions:
            samples = ",".join(failed_regions[region])
            failed_output.write("{}\t{}\n".format(region, samples))

    with open(args.depth_output, 'w') as depth_output:
        depth_output.write("Region")
        for sample in samples_list:
            depth_output.write("\t{}".format(sample))
        depth_output.write("\n")

        for region in depth_regions:
            depth_output.write("{}\t".format(region))
            for sample in samples_list:
                depth_output.write("\t{}".format(depth_regions[region][sample]))
            depth_output.write("\n")
