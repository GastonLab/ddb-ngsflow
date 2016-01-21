#!/usr/bin/env python

import os
import csv
import sys
import argparse
import pybedtools
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plot

from collections import defaultdict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="CSV file containing sample meta-data")
    parser.add_argument('-d', '--depth_output', help="Output file name for depth information")
    parser.add_argument('-f', '--failed_output', help="Output file name for regions with failed samples")
    parser.add_argument('-o', '--output_summary', help="Output containing summary data and statistics")
    args = parser.parse_args()

    failed_regions = defaultdict(list)
    samples_metadata = defaultdict(dict)
    regions_summary_stats = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: -1)))
    depth_regions = defaultdict(lambda: defaultdict(lambda: -1))
    samples_list = list()

    sys.stdout.write("Reading sample metadata\n")
    if args.samples_file:
        with open(args.samples_file, 'r') as samples_file:
            reader = csv.reader(samples_file, delimiter=',')
            reader.next()
            for row in reader:
                samples_metadata[row[0]]['Pool'] = row[4]
                samples_metadata[row[0]]['NumSamples'] = row[1]
                samples_metadata[row[0]]['Combo'] = row[2]
                samples_metadata[row[0]]['Extraction'] = row[3]

    sys.stdout.write("Reading coverage data files\n")
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

    sys.stdout.write("Parsing data by categories\n")
    data_lists = list()
    data_all = list()
    four_samples_lists = list()
    eight_samples_lists = list()
    sixteen_samples_lists = list()
    four_samples_all = list()
    eight_samples_all = list()
    sixteen_samples_all = list()
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
            if depth_regions[region][sample] is not -1:
                temp.append(float(depth_regions[region][sample]))

                # sys.stdout.write("DEBUG: Num Samples: {}\n".format(samples_metadata[sample]['NumSamples']))
                # sys.stdout.write("DEBUG: Sample: {}\n".format(sample))
                if int(samples_metadata[sample]['NumSamples']) == 4:
                    four_samples_temp.append(float(depth_regions[region][sample]))
                if int(samples_metadata[sample]['NumSamples']) == 8:
                    eight_samples_temp.append(float(depth_regions[region][sample]))
                if int(samples_metadata[sample]['NumSamples']) == 16:
                    sixteen_samples_temp.append(float(depth_regions[region][sample]))

        # sys.stdout.write("DEBUG: Full Region {}: {}\n".format(region, len(temp)))
        # sys.stdout.write("DEBUG: Four Samples Region {}: {}\n".format(region, len(four_samples_temp)))
        # sys.stdout.write("DEBUG: Eight Samples Region {}: {}\n".format(region, len(eight_samples_temp)))
        # sys.stdout.write("DEBUG: Sixteen Samples Region {}: {}\n".format(region, len(sixteen_samples_temp)))

        regions_summary_stats[region]['all']['mean'] = np.mean(temp)
        regions_summary_stats[region]['all']['median'] = np.median(temp)
        regions_summary_stats[region]['all']['std'] = np.std(temp)

        # regions_summary_stats[region]['four_samples']['mean'] = np.mean(four_samples_temp)
        # regions_summary_stats[region]['four_samples']['median'] = np.median(four_samples_temp)
        # regions_summary_stats[region]['four_samples']['std'] = np.std(four_samples_temp)
        #
        # regions_summary_stats[region]['eight_samples']['mean'] = np.mean(eight_samples_temp)
        # regions_summary_stats[region]['eight_samples']['median'] = np.median(eight_samples_temp)
        # regions_summary_stats[region]['eight_samples']['std'] = np.std(eight_samples_temp)
        #
        # regions_summary_stats[region]['sixteen_samples']['mean'] = np.mean(sixteen_samples_temp)
        # regions_summary_stats[region]['sixteen_samples']['median'] = np.median(sixteen_samples_temp)
        # regions_summary_stats[region]['sixteen_samples']['std'] = np.std(sixteen_samples_temp)

        data_lists.append(temp)
        four_samples_lists.append(four_samples_temp)
        eight_samples_lists.append(eight_samples_temp)
        sixteen_samples_lists.append(sixteen_samples_temp)

        data_all.extend(temp)
        four_samples_all.extend(four_samples_temp)
        eight_samples_all.extend(eight_samples_temp)
        sixteen_samples_all.extend(sixteen_samples_temp)

    data = np.array(data_lists)
    four_samples_data = np.array(four_samples_lists)
    eight_samples_data = np.array(eight_samples_lists)
    sixteen_samples_data = np.array(sixteen_samples_lists)

    regions_summary_stats['Total']['all']['mean'] = np.mean(data_all)
    regions_summary_stats['Total']['all']['median'] = np.median(data_all)
    regions_summary_stats['Total']['all']['std'] = np.std(data_all)

    # regions_summary_stats['Total']['four_samples']['mean'] = np.mean(four_samples_all)
    # regions_summary_stats['Total']['four_samples']['median'] = np.median(four_samples_all)
    # regions_summary_stats['Total']['four_samples']['std'] = np.std(four_samples_all)
    #
    # regions_summary_stats['Total']['eight_samples']['mean'] = np.mean(eight_samples_all)
    # regions_summary_stats['Total']['eight_samples']['median'] = np.median(eight_samples_all)
    # regions_summary_stats['Total']['eight_samples']['std'] = np.std(eight_samples_all)
    #
    # regions_summary_stats['Total']['sixteen_samples']['mean'] = np.mean(sixteen_samples_all)
    # regions_summary_stats['Total']['sixteen_samples']['median'] = np.median(sixteen_samples_all)
    # regions_summary_stats['Total']['sixteen_samples']['std'] = np.std(sixteen_samples_all)

    with open(args.output_summary, 'w') as summary:
        summary.write("Region\tAll Mean\tAll Median\tAll Std\tFour Mean\tFour Median\tFour Std\tEight Mean\t"
                      "Eight Median\t Eight Std\tSixteen Mean\tSixteen Median\tSixteen Std\n")
        for region in depth_regions:
            summary.write("{region}\t{allm}\t{allme}\t{alls}\t{fourm}\t{fourme}\t{fours}\t{eightm}\t{eightme}\t{eights}"
                          "\t{sixteenm}\t{sixteenme}\t"
                          "{sixteens}\n".format(region=region,
                                                allm=regions_summary_stats[region]['all']['mean'],
                                                allme=regions_summary_stats[region]['all']['median'],
                                                alls=regions_summary_stats[region]['all']['std'],
                                                fourm=regions_summary_stats[region]['four_samples']['mean'],
                                                fourme=regions_summary_stats[region]['four_samples']['median'],
                                                fours=regions_summary_stats[region]['four_samples']['std'],
                                                eightm=regions_summary_stats[region]['eight_samples']['mean'],
                                                eightme=regions_summary_stats[region]['eight_samples']['median'],
                                                eights=regions_summary_stats[region]['eight_samples']['std'],
                                                sixteenm=regions_summary_stats[region]['sixteen_samples']['mean'],
                                                sixteenme=regions_summary_stats[region]['sixteen_samples']['median'],
                                                sixteens=regions_summary_stats[region]['sixteen_samples']['std']))
        summary.write("Summary\t{allm}\t{allme}\t{alls}\t{fourm}\t{fourme}\t{fours}\t{eightm}\t{eightme}\t{eights}"
                      "\t{sixteenm}\t{sixteenme}\t"
                      "{sixteens}\n".format(allm=regions_summary_stats['Total']['all']['mean'],
                                            allme=regions_summary_stats['Total']['all']['median'],
                                            alls=regions_summary_stats['Total']['all']['std'],
                                            fourm=regions_summary_stats['Total']['four_samples']['mean'],
                                            fourme=regions_summary_stats['Total']['four_samples']['median'],
                                            fours=regions_summary_stats['Total']['four_samples']['std'],
                                            eightm=regions_summary_stats['Total']['eight_samples']['mean'],
                                            eightme=regions_summary_stats['Total']['eight_samples']['median'],
                                            eights=regions_summary_stats['Total']['eight_samples']['std'],
                                            sixteenm=regions_summary_stats['Total']['sixteen_samples']['mean'],
                                            sixteenme=regions_summary_stats['Total']['sixteen_samples']['median'],
                                            sixteens=regions_summary_stats['Total']['sixteen_samples']['std']))

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

    # count = 1
    # locations = list()
    # for name in names:
    #     locations.append(count)
    #     count += 1

    # plot.figure(1)
    # plot1 = plot.boxplot(data.transpose(), widths=0.7, positions=locations, notch=True, patch_artist=True)
    #
    # plot.figure(2)
    # plot2 = plot.boxplot(four_samples_data.transpose(), widths=0.7, positions=locations, notch=True, patch_artist=True)
    #
    # plot.figure(3)
    # plot3 = plot.boxplot(eight_samples_data.transpose(), widths=0.7, positions=locations, notch=True, patch_artist=True)
    #
    # plot.figure(4)
    # plot4 = plot.boxplot(sixteen_samples_data.transpose(), widths=0.7, positions=locations, notch=True, patch_artist=True)
    #
    # plot.xticks(locations, names, rotation='vertical')
    # plot.ylabel('Depth of Coverage')
    #
    # plot.show()
