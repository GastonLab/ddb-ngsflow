"""
.. module:: utilities
   :platform: Unix, OSX
   :synopsis: A module of methods for various utility functions. This is a general catch all module for methoods
   not better implemented elsewhere.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>

"""

import sys
import pyexcel
import pyexcel.ext.xlsx
import pybedtools
import multiprocessing

from collections import defaultdict

from ddb_ngsflow import pipeline


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


def sambamba_region_coverage(job, config, sample, samples, input_bam):
    """Run SamBambam to calculate the coverage of targeted regions
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type samples: dict
    :param samples: The samples configuration dictionary
    :type input_bam: str.
    :returns:  str -- The output BED file name.
    """

    output = "{}.sambamba_coverage.bed".format(sample)
    logfile = "{}.sambamba_coverage.log".format(sample)

    command = ("{}".format(config['sambamba']['bin']),
               "depth region",
               "-L",
               "{}".format(samples[sample]['region']),
               "-t",
               "{}".format(config['sambamba']['num_cores']),
               "-T".format(config['coverage_threshold']),
               "{}".format(input_bam),
               ">",
               "{}".format(output))

    job.fileStore.logToMaster("SamBamba Coverage Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return output


def bedtools_coverage_per_site(job, config, sample, input_bam):
    """Run BedTools to calculate the per-site coverage of targeted regions
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output BED file name.
    """

    output = "{}.bedtools_coverage_per_site.bed".format(sample)
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
    :type config: dict..
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


def read_coverage(job, config, sample, vcf):
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


def bcftools_filter_variants_regions(job, config, sample, samples, input_vcf):
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
                      "{}".format(samples[sample]['regions']),
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


def merge_samples(job, config, sample, input_vcf1, input_vcf2):
    """Merge samples into a single VCF"""

    output_vcf = "{}.ds.merged.vcf".format(sample)
    logfile = "{}.ds_merging.log".format(sample)

    bgzipped_vcf1 = "{}.gz".format(input_vcf1)
    bgzipped_vcf2 = "{}.gz".format(input_vcf2)

    bgzip_and_tabix_vcf(job, input_vcf1)
    bgzip_and_tabix_vcf(job, input_vcf2)

    # Prefer to use BCFTools but merge returns errors about LSEQ fields
    # command = ("{}".format(config['bcftools']['bin']),
    #            "merge",
    #            "-O",
    #            "v",
    #            "-o",
    #            "{}".format(output_vcf),
    #            "{}".format(bgzipped_vcf1),
    #            "{}".format(bgzipped_vcf2))

    command = ("{}".format(config['vcftools_merge']['bin']),
               "{}".format(bgzipped_vcf1),
               "{}".format(bgzipped_vcf2),
               ">",
               "{}".format(output_vcf))

    job.fileStore.logToMaster("Merging stranded libraries into multi-sample VCF with the command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return output_vcf


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
