"""
.. module:: utilities
   :platform: Unix, OSX
   :synopsis: A module of methods for various utility functions. This is a general catch all module for methoods
   not better implemented elsewhere.

.. moduleauthor:: Daniel Gaston <daniel.gaston@gmail.com>


"""

__author__ = 'dgaston'

import sys
import pyexcel
import pyexcel.ext.xlsx
import pybedtools
import multiprocessing

from ngsflow import pipeline


def spawn_batch_jobs(job):
    """
    This is simply a placeholder root job for the workflow
    """

    job.fileStore.logToMaster("Initializing workflow\n")


def spawn_variant_jobs(job):
    """
    This is simply a placeholder job to create a node in the graph for spawning
    off the multiple variant callers
    """

    job.fileStore.logToMaster("Spawning all variant calling methods\n")


def run_fastqc(job, config, samples):
    """Run FastQC on provided FastQ files
    :param config: The configuration dictionary.
    :type config: dict.
    :param samples: Samples dictionary
    :type samples: str.
    """

    job.fileStore.logToMaster("Running FastQC for all samples\n")
    logfile = "fastqc.log"

    fastq_files_list = list()
    for sample in samples:
        fastq_files_list.append(samples[sample]['fastq1'])
        fastq_files_list.append(samples[sample]['fastq2'])

    if multiprocessing.cpu_count() <= len(samples):
        num_cores = multiprocessing.cpu_count()
    else:
        num_cores = len(samples)
    fastq_files_string = " ".join(fastq_files_list)
    command = ("{}".format(config['fastqc']['bin']),
               "{}".format(fastq_files_string),
               "--extract",
               "-t",
               "{}".format(num_cores))

    job.fileStore.logToMaster("FastQC Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)


def generate_fastqc_summary_report(job, config, samples):
    """Parse FastQC summary reports and generate a run-level summary
    :param config: The configuration dictionary.
    :type config: dict.
    :param samples: Samples dictionary
    :type samples: str.
    """

    job.fileStore.logToMaster("Parsing FastQC results to run-level summary file\n")
    with open("{}_fastqc_summary.txt".format(config['run_name']), 'w') as summary_file:
        for sample in samples:
            sample_fastq_dirs = list()
            if sample['fastq1']:
                temp = sample['fastq1'].split(".")
                sample_fastq_dirs.append(temp[0])
            if sample['fastq2']:
                temp = sample['fastq2'].split(".")
                sample_fastq_dirs.append(temp[0])
            for dirbase in sample_fastq_dirs:
                with open("./%s_fastqc/summary.txt" % dirbase, "rU") as fastqc_file:
                    for line in fastqc_file.read():
                        summary_file.write(line)


def vt_normalization(job, config, sample, input_vcf):
    """Decompose and left normalize variants
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_vcf: The input_vcf file name to process.
    :type input_vcf: str.
    :returns:  str -- The output vcf file name.
    """

    output_vcf = "{}.normalized.vcf".format(sample)
    logfile = "{}.vt_normalization.log".format(sample)

    normalization = ("zless",
                     "{}".format(input_vcf),
                     "|",
                     "sed",
                     "'s/ID=AD,Number=./ID=AD,Number=R/'",
                     "|",
                     "vt",
                     "decompose",
                     "-s",
                     "-",
                     "|",
                     "vt",
                     "normalize",
                     "-r",
                     "{}".format(config['reference']),
                     "-",
                     ">",
                     "{}".format(output_vcf))

    job.fileStore.logToMaster("VT Command: {}\n".format(normalization))
    pipeline.run_and_log_command(" ".join(normalization), logfile)

    return output_vcf


def bedtools_coverage_per_site(job, config, sample, input_bam):
    """Run BedTools to calculate the per-site coverage of targeted regions
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_vcf file name to process.
    :type input_bam: str.
    :returns:  str -- The output BED file name.
    """

    output = "{}.coverage.bed".format(sample)
    logfile = "{}.bedtools_coverage.log".format(sample)

    coverage = ("{}".format(config['bedtools']['bin']),
                "coverage",
                "-d",
                "-a",
                "{}".format(config['regions']),
                "-b",
                "{}".format(input_bam),
                ">",
                "{}".format(output))

    job.fileStore.logToMaster("BedTools Coverage Command: {}\n".format(coverage))
    pipeline.run_and_log_command(" ".join(coverage), logfile)

    return output


def bedtools_coverage_to_summary(job, config, sample, input_file):
    """Summarize outputs from BedTools coverage results"""

    raise NotImplementedError


def generate_coverage_report(job, config, vcfs):
    """Take DiagnoseTargets data and generate a coverage report
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param vcfs: The list of input DiagnoseTargets generated vcf file names to process.
    :type vcfs: str.
    """

    samples_coverage = {"Chr": [], "Start": [], "Stop": [], "Target": []}
    first_pass = True

    job.fileStore.logToMaster("Processing DiagnoseTargets outputs and writing to spreadsheet\n")
    sys.stdout.write("Processing VCFs:\n")
    for vcf in vcfs:
        sys.stdout.write("{}\n".format(vcf))

    for vcf in vcfs:
        filter_field = "{}_filter".format(vcf)
        depth_field = "{}_depth".format(vcf)
        low_field = "{}_bp_low".format(vcf)
        zero_field = "{}_bp_zero".format(vcf)

        samples_coverage[filter_field] = list()
        samples_coverage[depth_field] = list()
        samples_coverage[low_field] = list()
        samples_coverage[zero_field] = list()

        targeted_regions = pybedtools.BedTool(config['regions'])
        coverage_data = pybedtools.BedTool(vcf)
        intersections = coverage_data.intersect(targeted_regions, loj=True)

        for region in intersections:
            if first_pass:
                samples_coverage['Chr'].append(region.chrom)
                samples_coverage['Start'].append(region.start)
                samples_coverage['Stop'].append(region.stop)
                samples_coverage['Target'].append(region[13])
            reads_data = region[9].split(":")
            samples_coverage[filter_field].append(region[6])
            samples_coverage[depth_field].append(reads_data[-3])
            samples_coverage[low_field].append(reads_data[-2])
            samples_coverage[zero_field].append(reads_data[-1])

        if first_pass:
            first_pass = False

    content = pyexcel.utils.dict_to_array(samples_coverage)
    sheet = pyexcel.Sheet(content)
    sheet.save_as("{}_coverage_results.xlsx".format(config['run_name']))


def bcftools_filter_variants_regions(job, config, sample, input_vcf):
    """Use bcftools to filter vcf file to only variants found within the specified regions file
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_vcf: The input_vcf file name to process.
    :type input_vcf: str.
    :returns:  str -- The output vcf file name.
    """

    filtered_vcf = "{}.on_target.vcf".format(sample)
    sorted_vcf = "{}.on_target_sorted.vcf".format(sample)
    bgzipped_vcf = "{}.gz".format(input_vcf)
    logfile = "{}.on_target_filter.log".format(sample)
    sort_logfile = "{}.on_target_sorted.log".format(sample)

    bgzip_and_tabix_vcf(job, input_vcf)

    filter_command = ("{}".format(config['bcftools']['bin']),
                      "isec",
                      "-T",
                      "{}".format(config['regions']),
                      "{}".format(bgzipped_vcf),
                      ">",
                      "{}".format(filtered_vcf))

    sort_command = ("cat",
                    "{}".format(filtered_vcf),
                    "|",
                    "{}".format(config['vcftools_sort']['bin']),
                    "-c",
                    ">",
                    "{}".format(sorted_vcf)
                    )

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
    """Run GATK's VariantFilter on the specified VCF
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
