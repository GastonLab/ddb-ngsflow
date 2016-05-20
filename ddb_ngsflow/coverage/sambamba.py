"""
.. module:: utilities
   :platform: Unix, OSX
   :synopsis: A module of methods for various sambamba functions. This is a general catch all module for methoods
   not better implemented elsewhere.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>

"""

from collections import defaultdict
from ddb_ngsflow import pipeline


def sambamba_region_coverage(job, config, name, samples, input_bam):
    """Run SamBambam to calculate the coverage of targeted regions
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample/library name.
    :type name: str.
    :param input_bam: The input_bam file name to process.
    :type samples: dict
    :param samples: The samples configuration dictionary
    :type input_bam: str.
    :returns:  str -- The output BED file name.
    """

    output = "{}.sambamba_coverage.bed".format(name)
    logfile = "{}.sambamba_coverage.log".format(name)

    command = ("{}".format(config['sambamba']['bin']),
               "depth region",
               "-L",
               "{}".format(samples[name]['regions']),
               "-t",
               "{}".format(config['sambamba']['num_cores']),
               "-T",
               "{}".format(config['coverage_threshold']),
               "-T",
               "{}".format(config['coverage_threshold2']),
               "{}".format(input_bam),
               ">",
               "{}".format(output))

    job.fileStore.logToMaster("SamBamba Coverage Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return output


def sambamba_coverage_analysis_job(job, config, samples, outfile):
    amplicon_coverage = defaultdict(lambda: defaultdict(int))
    num_samples = 0.0
    num_t = 0.0
    num_mp = 0.0
    num_cmp = 0.0
    for sample in samples:
        num_samples += 1
        if samples[sample]['extraction'] == 'T':
            num_t += 1
        elif samples[sample]['extraction'] == 'MP':
            num_mp += 1
        elif samples[sample]['extraction'] == 'CMP':
            num_cmp += 1
        else:
            sys.stderr.write("Unknown extraction protocol {}\n".format(samples[sample]['extraction']))
            sys.exit()
        coverage_file = "{}.sambamba_coverage.bed".format(sample)
        with open(coverage_file, 'rb') as coverage:
            reader = csv.reader(coverage, delimiter='\t')
            # Skip header
            reader.next()
            for row in reader:
                amplicon = "{}-{}-{}-{}".format(row[0], row[1], row[2], row[3])

                amplicon_coverage[amplicon][sample] = int(row[4])
                amplicon_coverage[amplicon]["{}_percent_{}".format(sample,
                                                                   config['coverage_threshold'])] = float(row[6])
                amplicon_coverage[amplicon]["{}_percent_{}".format(sample,
                                                                   config['coverage_threshold2'])] = float(row[7])
                amplicon_coverage[amplicon]['readcount_total'] += int(row[4])
                amplicon_coverage[amplicon]["readcount_{}_total".format(samples[sample]['extraction'])] += int(row[4])

                amplicon_coverage[amplicon]['percent_{}_total'.format(config['coverage_threshold'])] += float(row[6])

                amplicon_coverage[amplicon]['percent_{}_total'.format(config['coverage_threshold2'])] += float(row[7])
                amplicon_coverage[amplicon]['percent_{}_total_{}'.format(config['coverage_threshold'],
                                                                         samples[sample]['extraction'])] += float(row[6])
                amplicon_coverage[amplicon]['percent_{}_total_{}'.format(config['coverage_threshold2'],
                                                                         samples[sample]['extraction'])] += float(row[7])

                if float(row[6]) >= 100:
                    amplicon_coverage[amplicon]['num_samples_{}'.format(config['coverage_threshold'])] += 1
                    amplicon_coverage[amplicon]['num_samples_{}_{}'.format(config['coverage_threshold'],
                                                                           samples[sample]['extraction'])] += 1
                if float(row[7]) >= 100:
                    amplicon_coverage[amplicon]['num_samples_{}'.format(config['coverage_threshold2'])] += 1
                    amplicon_coverage[amplicon]['num_samples_{}_{}'.format(config['coverage_threshold2'],
                                                                           samples[sample]['extraction'])] += 1

    with open(outfile, 'wb') as summary:
        summary.write("Amplicon\tAvg Num Read per Sample\tPercent Samples {t1}\tPercent Samples {t2}\t"
                      "Percent Samples {t1} T\tPercent Samples {t2} T\t"
                      "Percent Samples {t1} MP\tPercent Samples {t2} MP\t"
                      "Percent Samples {t1} CMP\tPercent Samples {t2} CMP\t"
                      "Avg Percent bases at {t1}\tAvg Percent bases at {t2}\t"
                      "Avg Percent bases at {t1} T\tAvg Percent bases at {t2} T\t"
                      "Avg Percent bases at {t1} MP\tAvg Percent bases at {t2} MP\t"
                      "Avg Percent bases at {t1} CMP\tAvg Percent bases at {t2} CMP"
                      "\n".format(t1=config['coverage_threshold'],
                                  t2=config['coverage_threshold2']))

        for amplicon in amplicon_coverage:
            avg_reads = amplicon_coverage[amplicon]['readcount_total'] / num_samples

            perc_samples1 = amplicon_coverage[amplicon]['num_samples_{}'.format(config['coverage_threshold'])] / num_samples
            perc_samples2 = amplicon_coverage[amplicon]['num_samples_{}'.format(config['coverage_threshold2'])] / num_samples

            perc_samples_t1 = amplicon_coverage[amplicon]['num_samples_{}_{}'.format(config['coverage_threshold'],
                                                                                     "T")] / num_t
            perc_samples_t2 = amplicon_coverage[amplicon]['num_samples_{}_{}'.format(config['coverage_threshold2'],
                                                                                     "T")] / num_t

            perc_samples_mp1 = amplicon_coverage[amplicon]['num_samples_{}_{}'.format(config['coverage_threshold'],
                                                                                      "MP")] / num_mp
            perc_samples_mp2 = amplicon_coverage[amplicon]['num_samples_{}_{}'.format(config['coverage_threshold2'],
                                                                                      "MP")] / num_mp

            perc_samples_cmp1 = amplicon_coverage[amplicon]['num_samples_{}_{}'.format(config['coverage_threshold'],
                                                                                       "CMP")] / num_cmp
            perc_samples_cmp2 = amplicon_coverage[amplicon]['num_samples_{}_{}'.format(config['coverage_threshold2'],
                                                                                       "CMP")] / num_cmp

            avg_perc1 = amplicon_coverage[amplicon]['percent_{}_total'.format(config['coverage_threshold'])] / num_samples
            avg_perc2 = amplicon_coverage[amplicon]['percent_{}_total'.format(config['coverage_threshold2'])] / num_samples

            avg_perc_t1 = amplicon_coverage[amplicon]['percent_{}_total_{}'.format(config['coverage_threshold'],
                                                                                   "T")] / num_t
            avg_perc_t2 = amplicon_coverage[amplicon]['percent_{}_total_{}'.format(config['coverage_threshold2'],
                                                                                   "T")] / num_t

            avg_perc_mp1 = amplicon_coverage[amplicon]['percent_{}_total_{}'.format(config['coverage_threshold'],
                                                                                    "MP")] / num_mp
            avg_perc_mp2 = amplicon_coverage[amplicon]['percent_{}_total_{}'.format(config['coverage_threshold2'],
                                                                                    "MP")] / num_mp

            avg_perc_cmp1 = amplicon_coverage[amplicon]['percent_{}_total_{}'.format(config['coverage_threshold'],
                                                                                     "CMP")] / num_cmp
            avg_perc_cmp2 = amplicon_coverage[amplicon]['percent_{}_total_{}'.format(config['coverage_threshold2'],
                                                                                     "CMP")] / num_cmp

            summary.write("{amp}\t{avg_reads}\t{percs1}\t{percs2}\t{percst1}\t{percst2}\t{percsmp1}\t{percsmp2}\t"
                          "{percscmp1}\t{percscmp2}\t{avg_perc1}\t{avg_perc2}\t{avg_perct1}\t{avg_perct2}\t"
                          "{avg_percmp1}\t{avg_percmp2}\t"
                          "{avg_perccmp1}\t{avg_perccmp2}"
                          "\n".format(amp=amplicon, avg_reads=avg_reads, percs1=perc_samples1, percs2=perc_samples2,
                                      percst1=perc_samples_t1, percst2=perc_samples_t2,
                                      percsmp1=perc_samples_mp1, percsmp2=perc_samples_mp2,
                                      percscmp1=perc_samples_cmp1, percscmp2=perc_samples_cmp2,
                                      avg_perc1=avg_perc1, avg_perc2=avg_perc2,
                                      avg_perct1=avg_perc_t1, avg_perct2=avg_perc_t2,
                                      avg_percmp1=avg_perc_mp1, avg_percmp2=avg_perc_mp2,
                                      avg_perccmp1=avg_perc_cmp1, avg_perccmp2=avg_perc_cmp2))


def sambamba_coverage_summary(job, config, samples, summary_outfile, outfile):
    amplicon_coverage = defaultdict(lambda: defaultdict(int))
    num_samples = 0.0
    for sample in samples:
        num_samples += 1
        coverage_file = "{}.sambamba_coverage.bed".format(sample)
        with open(coverage_file, 'rb') as coverage:
            reader = csv.reader(coverage, delimiter='\t')
            # Skip header
            reader.next()
            for row in reader:
                amplicon = "{}-{}-{}-{}".format(row[0], row[1], row[2], row[3])

                amplicon_coverage[amplicon][sample] = int(row[4])
                amplicon_coverage[amplicon]["{}_percent_{}".format(sample, config['coverage_threshold'])] = float(row[6])
                amplicon_coverage[amplicon]["{}_percent_{}".format(sample, config['coverage_threshold2'])] = float(row[7])
                amplicon_coverage[amplicon]['readcount_total'] += int(row[4])
                amplicon_coverage[amplicon]['percent_{}_total'.format(config['coverage_threshold'])] += float(row[6])
                amplicon_coverage[amplicon]['percent_{}_total'.format(config['coverage_threshold2'])] += float(row[7])

                if float(row[6]) >= 100:
                    amplicon_coverage[amplicon]['num_samples_{}'.format(config['coverage_threshold'])] += 1
                if float(row[7]) >= 100:
                    amplicon_coverage[amplicon]['num_samples_{}'.format(config['coverage_threshold2'])] += 1

    with open(summary_outfile, 'wb') as summary:
        with open(outfile, 'wb') as output:
            summary.write("Amplicon\tAvg Num Read per Sample\tPercent Samples {t1}\tPercent Samples {t2}\t"
                          "Avg Percent bases at {t1}\t"
                          "Avg Percent bases at {t2}\n".format(t1=config['coverage_threshold'],
                                                               t2=config['coverage_threshold2']))

            output.write("Sample\tAmplicon\tNum Reads\tPercent Bases at {t1}\t"
                         "Percent Bases at {t2}\n".format(t1=config['coverage_threshold'],
                                                          t2=config['coverage_threshold2']))

            for amplicon in amplicon_coverage:
                avg_reads = amplicon_coverage[amplicon]['readcount_total'] / num_samples
                perc_samples1 = amplicon_coverage[amplicon]['num_samples_{}'.format(config['coverage_threshold'])] / num_samples
                perc_samples2 = amplicon_coverage[amplicon]['num_samples_{}'.format(config['coverage_threshold2'])] / num_samples

                avg_perc1 = amplicon_coverage[amplicon]['percent_{}_total'.format(config['coverage_threshold'])] / num_samples
                avg_perc2 = amplicon_coverage[amplicon]['percent_{}_total'.format(config['coverage_threshold2'])] / num_samples

                summary.write("{amp}\t{avg_reads}\t{percs1}\t{percs2}\t"
                              "{avg_perc1}\t{avg_perc2}\n".format(amp=amplicon, avg_reads=avg_reads,
                                                                  percs1=perc_samples1, percs2=perc_samples2,
                                                                  avg_perc1=avg_perc1, avg_perc2=avg_perc2))

                for sample in samples:
                    output.write("{sample}\t{amp}\t{samp_reads}\t{s_perc1}\t"
                                 "{s_perc1}\n".format(sample=sample,
                                                      amp=amplicon,
                                                      samp_reads=amplicon_coverage[amplicon][sample],
                                                      s_perc1=amplicon_coverage[amplicon]["{}_percent_{}".format(sample, config['coverage_threshold'])],
                                                      s_perc2=amplicon_coverage[amplicon]["{}_percent_{}".format(sample, config['coverage_threshold2'])]))
