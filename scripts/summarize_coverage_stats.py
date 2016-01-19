#!/usr/bin/env python

import os
import csv
import sys
import argparse
import pybedtools
import numpy as np
import matplotlib.pyplot as plot

from collections import defaultdict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="CSV file containing sample meta-data")
    parser.add_argument('-d', '--depth_output', help="Output file name for depth information")
    parser.add_argument('-f', '--failed_output', help="Output file name for regions with failed samples")
    args = parser.parse_args()

    failed_regions = defaultdict(list)
    samples_metadata = defaultdict(dict)
    depth_regions = defaultdict(lambda: defaultdict(lambda: -1))
    samples_list = list()

    with open(args.samples_file, 'r') as samples_file:
        reader = csv.reader(samples_file, delimiter=',')
        reader.next()
        for row in reader:
            samples_metadata[row[0]]['RunID'] = row[1]
            samples_metadata[row[0]]['TargetPool'] = row[2]
            samples_metadata[row[0]]['NumSamples'] = row[3]
            samples_metadata[row[0]]['Combo'] = row[4]
            samples_metadata[row[0]]['Extraction'] = row[5]

    files = os.listdir(os.getcwd())
    for coverage_file in files:
        sections = coverage_file.partition(".")
        if sections[2] == 'diagnosetargets.vcf':
            samples_list.append(sections[0])
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

    data_lists = list()
    four_samples_lists = list()
    eight_samples_lists = list()
    sixteen_samples_lists = list()
    combo_lists = list()
    extraction_protocol_lists = list()
    names = list()

    for region in depth_regions:
        names.append(region)

        temp = list()
        four_samples_temp = list()
        eight_samples_temp = list()
        sixteen_samples_temp = list()
        combo_temp = list()
        extract_temp = list()

        for sample in samples_list:
            if depth_regions[region][sample] is not None:
                temp.append(float(depth_regions[region][sample]))
                # if samples_metadata[sample]['NumSamples'] == 4:
                #     four_samples_temp.append(float(depth_regions[region][sample]))
                # if samples_metadata[sample]['NumSamples'] == 8:
                #     four_samples_temp.append(float(depth_regions[region][sample]))
                # if samples_metadata[sample]['NumSamples'] == 16:
                #     four_samples_temp.append(float(depth_regions[region][sample]))

        # sys.stdout.write("DEBUG: Region {}: {}\n".format(region, len(temp)))
        data_lists.append(temp)
        # four_samples_lists.append(four_samples_temp)
        # eight_samples_lists.append(eight_samples_temp)
        # sixteen_samples_lists.append(sixteen_samples_temp)

    # count = 1
    # locations = list()
    # for name in names:
    #     locations.append(count)
    #     count += 1
    #
    # data = np.array(data_lists)
    # four_samples_data = np.array(four_samples_lists)
    # eight_samples_data = np.array(eight_samples_lists)
    # sixteen_samples_data = np.array(sixteen_samples_lists)
    #
    # plot.figure()
    #
    # plot1 = plot.boxplot(data.transpose(), widths=0.7, positions=locations, notch=True, patch_artist=True)
    # plot2 = plot.boxplot(four_samples_data.transpose(), widths=0.7, positions=locations, notch=True, patch_artist=True)
    # plot3 = plot.boxplot(eight_samples_data.transpose(), widths=0.7, positions=locations, notch=True, patch_artist=True)
    # plot4 = plot.boxplot(sixteen_samples_data.transpose(), widths=0.7, positions=locations, notch=True, patch_artist=True)
    #
    # plot.xticks(locations, names, rotation='vertical')
    # plot.ylabel('Depth of Coverage')
    #
    # plot.show()

    with open(args.failed_output, 'w') as failed_output:
        for region in failed_regions:
            samples = ",".join(failed_regions[region])
            failed_output.write("{}\t{}\t{}\n".format(region, len(failed_regions[region]), samples))

    with open(args.depth_output, 'w') as depth_output:
        depth_output.write("Sample")
        for region in depth_regions:
            depth_output.write("\t{}".format(region))
        depth_output.write("\n")
        for sample in samples_list:
            depth_output.write("{}".format(sample))
            for region in depth_regions:
                depth_output.write("\t{}".format(float(depth_regions[region][sample])))
            depth_output.write("\n")
