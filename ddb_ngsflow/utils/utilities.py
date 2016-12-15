"""
.. module:: utilities
   :platform: Unix, OSX
   :synopsis: A module of methods for various utility functions. This is a general catch all module for methoods
   not better implemented elsewhere.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>

"""

import csv
import sys
import pybedtools
import multiprocessing

from collections import defaultdict

from ddb_ngsflow import pipeline


def readsort_bam(job, config, name, samples):
    pass


def bedtools_coverage_per_site(job, config, name, input_bam):
    """Run BedTools to calculate the per-site coverage of targeted regions
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output BED file name.
    """

    output = "{}.bedtools_coverage_per_site.bed".format(name)
    logfile = "{}.bedtools_coverage.log".format(name)

    coverage = ["{}".format(config['bedtools']['bin']),
                "coverage",
                "-d",
                "-a",
                "{}".format(config['regions']),
                "-b",
                "{}".format(input_bam),
                ">",
                "{}".format(output)]

    job.fileStore.logToMaster("BedTools Coverage Command: {}\n".format(coverage))
    pipeline.run_and_log_command(" ".join(coverage), logfile)

    return output


def bedtools_coverage_to_summary(job, config, name, input_file):
    """Summarize outputs from BedTools coverage results"""

    raise NotImplementedError


def read_coverage(job, config, name, vcf):
    """Take DiagnoseTargets data return summarized results
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param vcf: The input DiagnoseTargets generated vcf file name to process.
    :type vcf: str.
    """

    sample_coverage = defaultdict(dict)

    targeted_regions = pybedtools.BedTool(config['regions'])
    coverage_data = pybedtools.BedTool(vcf)
    intersections = coverage_data.intersect(targeted_regions, loj=True)

    for region in intersections:
        sample_coverage[region[13]]['filter_field'] = region[6]
        reads_data = region[9].split(":")
        sample_coverage[region[13]]['filter_field'] = region[6]
        sample_coverage[region[13]]['depth_field'] = reads_data[-3]
        sample_coverage[region[13]]['low_field'] = reads_data[-2]
        sample_coverage[region[13]]['zero_field'] = reads_data[-1]

    return sample_coverage


def generate_coverage_summary(job, config, samples):
    """Take Summarized DiagnoseTargets data and generate a coverage summary
    :param config: The configuration dictionary.
    :type config: dict.
    :param samples: summarized sample results.
    :type samples: dict.
    """
    with open("sample_coverage_summary.txt", 'w') as outfile:
        for sample in samples:
            for target in samples[sample].keys():
                targets_list = target.split("\t")
                target_string = ",".join(targets_list)
                if 'COVERAGE' in samples[sample][target]['filter_field'] or 'NO_READS' in samples[sample][target]['filter_field']:
                    outfile.write("{sample}\t{region}\t{filter}\n".format(sample=sample, region=target_string, filter=samples[sample][target]['filter_field']))


def bcftools_filter_variants_regions(job, config, name, samples, input_vcf):
    """Use bcftools to filter vcf file to only variants found within the specified regions file
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_vcf: The input_vcf file name to process.
    :type input_vcf: str.
    :returns:  str -- The output vcf file name.
    """

    filtered_vcf = "{}.on_target.vcf".format(name)
    sorted_vcf = "{}.on_target_sorted.vcf".format(name)
    bgzipped_vcf = "{}.gz".format(input_vcf)
    logfile = "{}.on_target_filter.log".format(name)
    sort_logfile = "{}.on_target_sorted.log".format(name)

    bgzip_and_tabix_vcf(job, input_vcf)

    filter_command = ["{}".format(config['bcftools']['bin']),
                      "isec",
                      "-T",
                      "{}".format(samples[name]['regions']),
                      "{}".format(bgzipped_vcf),
                      ">",
                      "{}".format(filtered_vcf)]

    sort_command = ["cat",
                    "{}".format(filtered_vcf),
                    "|",
                    "{}".format(config['vcftools_sort']['bin']),
                    "-c",
                    ">",
                    "{}".format(sorted_vcf)
                    ]

    job.fileStore.logToMaster("BCFTools isec command for filtering to only target regions: {}\n".format(filter_command))
    pipeline.run_and_log_command(" ".join(filter_command), logfile)

    job.fileStore.logToMaster("VCFTools-sort command for filtering to only target regions: {}\n".format(sort_command))
    pipeline.run_and_log_command(" ".join(sort_command), sort_logfile)

    return sorted_vcf


def _bgzip_and_tabix_vcf_instructions(infile):
    """Generate instructions and logfile for bgzip and tabix"""

    bgzip_command = "bgzip -c %s > %s.gz" % (infile, infile)
    bgzip_logfile = "%s.bgzip.log" % infile

    tabix_command = "tabix -p vcf %s.gz" % infile
    tabix_logfile = "%s.tabix.log" % infile

    bgzip_instructions = list()
    bgzip_instructions.append(bgzip_command)
    bgzip_instructions.append(bgzip_logfile)

    tabix_instructions = list()
    tabix_instructions.append(tabix_command)
    tabix_instructions.append(tabix_logfile)

    return bgzip_instructions, tabix_instructions


def bgzip_and_tabix_vcf(job, infile):
    """Run BGZip and Tabix on the specified VCF
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param infile: The input_vcf file name to process.
    :type infile: str.
    :returns:  str -- The output vcf file name.
    """

    bgzip_instructions, tabix_instructions = _bgzip_and_tabix_vcf_instructions(infile)

    job.fileStore.logToMaster("BGzip Command: {}\n".format(bgzip_instructions[0]))
    pipeline.run_and_log_command(bgzip_instructions[0], bgzip_instructions[1])

    job.fileStore.logToMaster("Tabix Command: {}\n".format(tabix_instructions[0]))
    pipeline.run_and_log_command(tabix_instructions[0], tabix_instructions[1])
