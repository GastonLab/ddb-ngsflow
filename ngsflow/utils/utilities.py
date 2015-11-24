__author__ = 'dgaston'

import os
import sys
import pyexcel
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
    """Run FastQC on provided FastQ files"""

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
    """Parse FastQC summary reports and generate a run-level summary"""

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
    """Decompose and left normalize variants"""

    output_vcf = "{}.normalized.vcf".format(sample)
    logfile = "{}.vt_normalization.log".format(sample)

    normalization = ("zless",
                     "{}".format(input_vcf),
                     "|",
                     "sed",
                     "'s/ID=AD,Number=./ID=AD,Number=R/'",
                     "|",
                     "vt",
                     "decompose"
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
    """Run BedTools to calculate the per-site coverage of targeted regions"""

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


def generate_coverage_report(project_config, sample_config, tool_config, resource_config):
    """Take DiagnoseTargets data and generate a coverage report"""

    project_config['coverage_problems'] = "%s.coverage_issues_report.txt" % project_config['project_name']
    samples_coverage = {"Chr": [], "Start": [], "Stop": [], "Target": []}
    first_pass = True

    for sample in sample_config:
        filter_field = "%s_filter" % sample['name']
        depth_field = "%s_depth" % sample['name']
        low_field = "%s_bp_low" % sample['name']
        zero_field = "%s_bp_zero" % sample['name']

        samples_coverage[filter_field] = []
        samples_coverage[depth_field] = []
        samples_coverage[low_field] = []
        samples_coverage[zero_field] = []

        targeted_regions = pybedtools.BedTool(resource_config['regions'])
        coverage_data = pybedtools.BedTool(sample['diagnose_targets_vcf'])
        intersections = coverage_data.intersect(targeted_regions, loj=True)

        for region in intersections:
            if first_pass:
                samples_coverage['Chr'].append(region.chrom)
                samples_coverage['Start'].append(region.start)
                samples_coverage['Stop'].append(region.stop)
                samples_coverage['Target'].append(region[13])
            reads_data = region[9].split(":")
            samples_coverage[filter_field].append(region[6])
            samples_coverage[depth_field].append(reads_data[0])
            samples_coverage[low_field].append(reads_data[1])
            samples_coverage[zero_field].append(reads_data[2])

        if first_pass:
            first_pass = False

    content = pyexcel.utils.dict_to_array(samples_coverage)
    sheet = pyexcel.Sheet(content)
    sheet.save_as("%s_coverage_results.xlsx" % project_config['project_name'])


def bcftools_filter_variants_regions(job, config, sample, input_vcf):
    """Use bcftools to filter vcf file to only variants found within the specified regions file"""

    filtered_vcf = "{}.on_target.vcf".format(sample)
    bgzipped_vcf = "{}.gz".format(input_vcf)
    logfile = "{}.on_target_filter.log".format(sample)

    bgzip_and_tabix_vcf(job, input_vcf)

    filter_command = ("{}".format(config['bcftools']['bin']),
                      "isec",
                      "-T",
                      "{}".format(config['regions']),
                      "{}".format(bgzipped_vcf),
                      ">",
                      "{}".format(filtered_vcf))

    job.fileStore.logToMaster("BCFTools isec command for filtering to only target regions: {}\n".format(filter_command))
    pipeline.run_and_log_command(" ".join(filter_command), logfile)

    return filtered_vcf


def bgzip_and_tabix_vcf_instructions(infile):
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
    """Call bgzip and tabix on vcf files"""

    bgzip_instructions, tabix_instructions = bgzip_and_tabix_vcf_instructions(infile)

    job.fileStore.logToMaster("BGzip Command: {}\n".format(bgzip_instructions[0]))
    pipeline.run_and_log_command(bgzip_instructions[0], bgzip_instructions[1])

    job.fileStore.logToMaster("Tabix Command: {}\n".format(tabix_instructions[0]))
    pipeline.run_and_log_command(tabix_instructions[0], tabix_instructions[1])
