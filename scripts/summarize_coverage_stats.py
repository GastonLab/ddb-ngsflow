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
    parser.add_argument('-p', '--double_stranded', help="Flag for whether this is a double stranded, pooled experiment",
                        type=bool)
    args = parser.parse_args()

    # Define depth thresholds
    critical = 250
    low = 500
    warn = 1000

    failed_regions = defaultdict(list)
    samples_metadata = defaultdict(dict)
    regions_summary_stats = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: -1)))
    depth_regions = defaultdict(lambda: defaultdict(lambda: -1))
    sample_numbers = defaultdict(lambda: 0)
    samples_list = list()

    sys.stdout.write("Reading sample metadata\n")
    if args.samples_file:
        with open(args.samples_file, 'r') as samples_file:
            reader = csv.reader(samples_file, delimiter=',')
            reader.next()
            for row in reader:
                samples_list.append(row[0])
                samples_metadata[row[0]]['Pool'] = row[4]
                samples_metadata[row[0]]['NumSamples'] = row[1]
                samples_metadata[row[0]]['Combo'] = row[2]
                samples_metadata[row[0]]['Extraction'] = row[3]

                # Collect the total number of samples or libraries for a specific number of samples in a run
                # or for different extraction protocols
                sample_numbers[row[1]] += 1
                sample_numbers[row[3]] += 1

    sys.stdout.write("Reading coverage data files\n")
    files = os.listdir(os.getcwd())
    for coverage_file in files:
        sections = coverage_file.partition(".")
        if sections[2] in ('diagnosetargets.vcf', 'snv_regions.diagnosetargets.vcf',
                           'indel_regions.diagnosetargets.vcf'):
            coverage_data = pybedtools.BedTool(coverage_file)
            for region in coverage_data:
                info = region[7].split(';')
                data = region[9].split(":")
                filter_field = region[6]

                # Assumes the depth (IDP) is third position from the end. Should fix this to use the format
                # field as a map.
                end = info[0].split('=')
                depth = float(data[-3])

                if args.double_stranded:
                    data2 = region[10].split(":")
                    depth += float(data2[-3])

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

    mp_extract_lists = list()
    cmp_extract_lists = list()
    t_extract_lists = list()

    mp_extract_all = list()
    cmp_extract_all = list()
    t_extract_all = list()

    names = list()

    num_mp = 0
    num_mp_critical = 0
    num_mp_low = 0
    num_mp_warn = 0

    num_cmp = 0
    num_cmp_critical = 0
    num_cmp_low = 0
    num_cmp_warn = 0

    num_t = 0
    num_t_critical = 0
    num_t_low = 0
    num_t_warn = 0

    for region in depth_regions:
        names.append(region)

        temp = list()

        four_samples_temp = list()
        eight_samples_temp = list()
        sixteen_samples_temp = list()

        mp_extract_temp = list()
        cmp_extract_temp = list()
        t_extract_temp = list()

        for sample in samples_list:
            if depth_regions[region][sample] is not -1:
                temp.append(float(depth_regions[region][sample]))

                # sys.stdout.write("DEBUG: Num Samples: {}\n".format(samples_metadata[sample]['NumSamples']))
                # sys.stdout.write("DEBUG: Sample: {}\n".format(sample))
                if int(samples_metadata[sample]['NumSamples']) == 4:
                    four_samples_temp.append(float(depth_regions[region][sample]))

                    if samples_metadata[sample]['Extraction'] == 'MP':
                        mp_extract_temp.append(float(depth_regions[region][sample]))
                    if samples_metadata[sample]['Extraction'] == 'CMP':
                        cmp_extract_temp.append(float(depth_regions[region][sample]))
                    if samples_metadata[sample]['Extraction'] == 'T':
                        t_extract_temp.append(float(depth_regions[region][sample]))

                if int(samples_metadata[sample]['NumSamples']) == 8:
                    eight_samples_temp.append(float(depth_regions[region][sample]))

                    if samples_metadata[sample]['Extraction'] == 'MP':
                        mp_extract_temp.append(float(depth_regions[region][sample]))
                    if samples_metadata[sample]['Extraction'] == 'CMP':
                        cmp_extract_temp.append(float(depth_regions[region][sample]))
                    if samples_metadata[sample]['Extraction'] == 'T':
                        t_extract_temp.append(float(depth_regions[region][sample]))

                if int(samples_metadata[sample]['NumSamples']) == 16:
                    sixteen_samples_temp.append(float(depth_regions[region][sample]))

                if samples_metadata[sample]['Extraction'] == 'MP':
                    num_mp += 1
                    if float(depth_regions[region][sample]) <= critical:
                        num_mp_critical += 1
                    elif float(depth_regions[region][sample]) <= low:
                        num_mp_low += 1
                    elif float(depth_regions[region][sample]) <= warn:
                        num_mp_warn += 1
                if samples_metadata[sample]['Extraction'] == 'CMP':
                    num_cmp += 1
                    if float(depth_regions[region][sample]) <= critical:
                        num_cmp_critical += 1
                    elif float(depth_regions[region][sample]) <= low:
                        num_cmp_low += 1
                    elif float(depth_regions[region][sample]) <= warn:
                        num_cmp_warn += 1
                if samples_metadata[sample]['Extraction'] == 'T':
                    num_t += 1
                    if float(depth_regions[region][sample]) <= critical:
                        num_t_critical += 1
                    elif float(depth_regions[region][sample]) <= low:
                        num_t_low += 1
                    elif float(depth_regions[region][sample]) <= warn:
                        num_t_warn += 1

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

        regions_summary_stats[region]['mp_extraction']['mean'] = np.mean(mp_extract_temp)
        regions_summary_stats[region]['mp_extraction']['median'] = np.median(mp_extract_temp)
        regions_summary_stats[region]['mp_extraction']['std'] = np.std(mp_extract_temp)

        regions_summary_stats[region]['cmp_extraction']['mean'] = np.mean(cmp_extract_temp)
        regions_summary_stats[region]['cmp_extraction']['median'] = np.median(cmp_extract_temp)
        regions_summary_stats[region]['cmp_extraction']['std'] = np.std(cmp_extract_temp)

        regions_summary_stats[region]['t_extraction']['mean'] = np.mean(t_extract_temp)
        regions_summary_stats[region]['t_extraction']['median'] = np.median(t_extract_temp)
        regions_summary_stats[region]['t_extraction']['std'] = np.std(t_extract_temp)

        data_lists.append(temp)
        four_samples_lists.append(four_samples_temp)
        eight_samples_lists.append(eight_samples_temp)
        sixteen_samples_lists.append(sixteen_samples_temp)

        mp_extract_lists.append(mp_extract_temp)
        cmp_extract_lists.append(cmp_extract_temp)
        t_extract_lists.append(t_extract_temp)

        data_all.extend(temp)
        four_samples_all.extend(four_samples_temp)
        eight_samples_all.extend(eight_samples_temp)
        sixteen_samples_all.extend(sixteen_samples_temp)

        mp_extract_all.extend(mp_extract_temp)
        cmp_extract_all.extend(cmp_extract_temp)
        t_extract_all.extend(t_extract_temp)

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

    regions_summary_stats['Total']['mp_extraction']['mean'] = np.mean(mp_extract_all)
    regions_summary_stats['Total']['mp_extraction']['median'] = np.median(mp_extract_all)
    regions_summary_stats['Total']['mp_extraction']['std'] = np.std(mp_extract_all)

    regions_summary_stats['Total']['cmp_extraction']['mean'] = np.mean(cmp_extract_all)
    regions_summary_stats['Total']['cmp_extraction']['median'] = np.median(cmp_extract_all)
    regions_summary_stats['Total']['cmp_extraction']['std'] = np.std(cmp_extract_all)

    regions_summary_stats['Total']['t_extraction']['mean'] = np.mean(t_extract_all)
    regions_summary_stats['Total']['t_extraction']['median'] = np.median(t_extract_all)
    regions_summary_stats['Total']['t_extraction']['std'] = np.std(t_extract_all)

    with open(args.output_summary, 'w') as summary:
        summary.write("Region\tAll Mean\tAll Median\tAll Std\t"
                      "Four Mean\tFour Median\tFour Std\t"
                      "Eight Mean\tEight Median\t Eight Std\t"
                      "Sixteen Mean\tSixteen Median\tSixteen Std\t"
                      "MP Mean\tMP Median\tMP Std\t"
                      "CMP Mean\tCMP Median\tCMP Std\t"
                      "T Mean\tT Median\tT Std\n")
        for region in depth_regions:

            summary.write("{region}\t{allm}\t{allme}\t{alls}\t"
                          "{fourm}\t{fourme}\t{fours}\t"
                          "{eightm}\t{eightme}\t{eights}\t"
                          "{sixteenm}\t{sixteenme}\t{sixteens}\t"
                          "{mpm}\t{mpme}\t{mps}\t"
                          "{cmpm}\t{cmpme}\t{cmps}\t"
                          "{tm}\t{tme}\t{ts}\t"
                          "\n".format(region=region,
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
                                      sixteens=regions_summary_stats[region]['sixteen_samples']['std'],
                                      mpm=regions_summary_stats[region]['mp_extraction']['mean'],
                                      mpme=regions_summary_stats[region]['mp_extraction']['median'],
                                      mps=regions_summary_stats[region]['mp_extraction']['std'],
                                      cmpm=regions_summary_stats[region]['cmp_extraction']['mean'],
                                      cmpme=regions_summary_stats[region]['cmp_extraction']['median'],
                                      cmps=regions_summary_stats[region]['cmp_extraction']['std'],
                                      tm=regions_summary_stats[region]['t_extraction']['mean'],
                                      tme=regions_summary_stats[region]['t_extraction']['median'],
                                      ts=regions_summary_stats[region]['t_extraction']['std']
                                      ))

        summary.write("Summary\t{allm}\t{allme}\t{alls}\t"
                      "{fourm}\t{fourme}\t{fours}\t"
                      "{eightm}\t{eightme}\t{eights}\t"
                      "{sixteenm}\t{sixteenme}\t{sixteens}\t"
                      "{mpm}\t{mpme}\t{mps}\t"
                      "{cmpm}\t{cmpme}\t{cmps}\t"
                      "{tm}\t{tme}\t{ts}\t"
                      "\n".format(allm=regions_summary_stats['Total']['all']['mean'],
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
                                  sixteens=regions_summary_stats['Total']['sixteen_samples']['std'],
                                  mpm=regions_summary_stats['Total']['mp_extraction']['mean'],
                                  mpme=regions_summary_stats['Total']['mp_extraction']['median'],
                                  mps=regions_summary_stats['Total']['mp_extraction']['std'],
                                  cmpm=regions_summary_stats['Total']['cmp_extraction']['mean'],
                                  cmpme=regions_summary_stats['Total']['cmp_extraction']['median'],
                                  cmps=regions_summary_stats['Total']['cmp_extraction']['std'],
                                  tm=regions_summary_stats['Total']['t_extraction']['mean'],
                                  tme=regions_summary_stats['Total']['t_extraction']['median'],
                                  ts=regions_summary_stats['Total']['t_extraction']['std']
                                  ))

    # with open(args.failed_output, 'w') as failed_output:
    #     for region in failed_regions:
    #         samples = ",".join(failed_regions[region])
    #         failed_output.write("{}\t{}\t{}\n".format(region, len(failed_regions[region]), samples))
    #
    # with open(args.depth_output, 'w') as depth_output:
    #     depth_output.write("Sample")
    #     for region in depth_regions:
    #         depth_output.write("\t{}".format(region))
    #     depth_output.write("\n")
    #     for sample in samples_list:
    #         depth_output.write("{}".format(sample))
    #         for region in depth_regions:
    #             depth_output.write("\t{}".format(float(depth_regions[region][sample])))
    #         depth_output.write("\n")

    with open("Extraction_Failure_Counts.txt", 'w') as extraction_stats:
        extraction_stats.write("Extraction Protocol (Total Number)\tNum Fail Critical ( <= {})\t"
                               "Num Fail Low ( <= {})\tNum Fail Warn ( <= {})\n".format(critical, low, warn))
        extraction_stats.write("MP ({})\t".format(num_mp))
        extraction_stats.write("{}\t{}\t{}\n".format(num_mp_critical, num_mp_low, num_mp_warn))

        extraction_stats.write("CMP ({})\t".format(num_cmp))
        extraction_stats.write("{}\t{}\t{}\n".format(num_cmp_critical, num_cmp_low, num_cmp_warn))

        extraction_stats.write("T ({})\t".format(num_t))
        extraction_stats.write("{}\t{}\t{}\n".format(num_t_critical, num_t_low, num_t_warn))
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
