__author__ = 'dgaston'

# This file contains methods for setting up and calling GATK and Picard tool executables in both per_sample and
# per_cohort modes. Haplotype caller, UnifiedGenotyper and MuTect are located in the variation.py file

import sys
import time
import multiprocessing

from ngsflow.utils import utilities
from ngsflow import pipeline


def run_diagnosetargets(project_config, sample_config, tool_config, resource_config):
    """Run GATK's DiagnoseTargets against the supplied regions"""

    instructions = list()
    command_core = ("java -Xmx%sg -jar %s -T DiagnoseTargets -R %s -L %s --minimum_coverage %s "
                    "--coverage_status_threshold 0.001" %
                    (tool_config['gatk']['max_mem'], tool_config['gatk']['bin'],
                     resource_config['reference_genome'], resource_config['regions'],
                     resource_config['coverage_threshold']))

    for sample in sample_config:
        sample['diagnose_targets_vcf'] = "%s.diagnosetargets.vcf" % sample['name']
        sample['diagnose_targets_missing_intervals'] = "%s.missing.intervals" % sample['name']
        logfile = "%s.diagnosetargets.log" % sample['name']

        command = ("%s -I %s -o %s --missing_intervals %s" % (command_core, sample['working_bam'],
                                                              sample['diagnose_targets_vcf'],
                                                              sample['diagnose_targets_missing_intervals']))
        instructions.append((command, logfile))

    sys.stdout.write("Running DiagnoseTargets\n")
    pipe.execute_multiprocess(instructions, int(tool_config['gatk']['num_cores']))
    sys.stdout.write("Running DiagnoseTargets\n")


def run_qualifymissing(project_config, sample_config, tool_config, resource_config):
    """Run GATK's QualifyMissingIntervals against the supplied regions and datae from DiagnoseTargets"""

    instructions = list()
    command_core = ("java -Xmx%sg -jar %s -T QualifyMissingIntervals -R %s" % (tool_config['gatk']['max_mem'],
                                                                               tool_config['gatk']['bin'],
                                                                               resource_config['reference_genome']))

    for sample in sample_config:
        sample['diagnose_targets_vcf'] = "%s.diagnosetargets.vcf" % sample['name']
        logfile = "%s.diagnosetargets.log" % sample['name']

        command = ("%s -I %s -L %s -o %s" % (command_core, sample['working_bam'],
                                             sample['diagnose_targets_missing_intervals'],
                                             sample['diagnose_targets_vcf']))
        instructions.append((command, logfile))

    sys.stdout.write("Running DiagnoseTargets\n")
    pipe.execute_multiprocess(instructions, int(tool_config['gatk']['num_cores']))
    sys.stdout.write("Running DiagnoseTargets\n")


def annotate_vcf(job, config, sample, input_vcf, input_bam, max_mem):
    """GATK Annotate and Variant Filters"""

    output_vcf = "{}.annotated.vcf".format(sample)
    annotation_log = "{}.variantannotation.log".format(sample)

    annotation_command = ("java",
                          "-Xmx{}g".format(max_mem),
                          "-jar",
                          "{}".format(config['gatk']),
                          "-T",
                          "VariantAnnotator",
                          "-R",
                          "{}".format(config['reference']),
                          "-nt",
                          "{}".format(multiprocessing.cpu_count()),
                          "--group",
                          "StandardAnnotation",
                          "--dbsnp"
                          "{}".format(config['dbsnp']),
                          "-I",
                          "{}".format(input_bam),
                          "--variant",
                          "{}".format(input_vcf),
                          "-L",
                          "{}".format(input_vcf),
                          "-o",
                          "{}".format(output_vcf))

    job.fileStore.logToMaster("GATK VariantAnnotator Command: {}\n".format(annotation_command))
    pipeline.run_and_log_command(" ".join(annotation_command), annotation_logfile)

    return output_vcf


def filter_variants(job, config, sample, input_vcf, max_mem):
    """Run GATK Variant Filtration"""

    output_vcf = "{}.filtered.vcf".format(sample)
    filter_log = "{}.variantfiltration.log".format(sample)

    filter_command = ("java",
                      "-Xmx{}g".format(max_mem),
                      "-jar",
                      "{}".format(config['gatk']),
                      "-T",
                      "VariantAnnotator",
                      "-R",
                      "{}".format(config['reference']),
                      "--filterExpression",
                      "'MQ0 > 50'",
                      "--filterName",
                      "'HighMQ0'",
                      "--filterExpression",
                      "'DP < {}'".format(config['coverage_threshold']),
                      "--filterName",
                      "'LowDepth'",
                      "--filterExpression",
                      "'QUAL < 10'",
                      "--filterName",
                      "'LowQual'",
                      "--filterExpression",
                      "'MQ < 10'",
                      "--filterName",
                      "'LowMappingQual'",
                      "--variant",
                      "{}".format(input_vcf),
                      "-o",
                      "{}".format(output_vcf))

    job.fileStore.logToMaster("GATK VariantFiltration Command: {}\n".format(filter_command))
    pipeline.run_and_log_command(" ".join(filter_command), filter_log)

    return output_vcf


def run_mark_duplicates(project_config, sample_config, tool_config, resource_config):
    """Run Picard MarkDuplicates"""

    sys.stdout.write("Running MarkDuplicates\n")
    instructions = list()
    for sample in sample_config:
        logfile = "%s.markduplicates.log" % sample['name']

        try:
            sample['sorted_bam'] = sample['bam']
        except KeyError:
            pass

        sample['dedup_bam'] = "%s.dedup.sorted.bam" % sample['name']
        metrics = "%s.dedup.metrics" % sample['name']

        command = ("java -Xmx%sg -jar %s MarkDuplicates CREATE_INDEX=true INPUT=%s OUTPUT=%s METRICS_FILE=%s "
                   "VALIDATION_STRINGENCY=LENIENT" %
                   (tool_config['gatk']['max_mem'], tool_config['picard']['bin'], sample['working_bam'],
                    sample['dedup_bam'], metrics))
        instructions.append((command, logfile))
        sample['working_bam'] = sample['dedup_bam']


def add_or_replace_readgroups(job, config, sample, input_bam, max_mem):
    """Run AddOrReplaceReadGroups"""

    job.fileStore.logToMaster("Running AddOrReplaceReadGroups in sample: {}".format(sample))

    output_bam = "{}.rg.sorted.bam".format(sample)
    logfile = "{}.addreadgroups.log".format(sample)
    index_log = "{}.buildindex.log".format(sample)

    command = ("java",
               "-Xmx{}g".format(max_mem),
               "-jar",
               "{}".format(config['picard']),
               "AddOrReplaceReadGroups",
               "INPUT={}".format(input_bam),
               "OUTPUT={}".format(output_bam),
               "RGID={}".format(sample),
               "RGSM={}".format(sample),
               "RGLB={}".format(sample),
               "RGPL=illumina",
               "RGPU=miseq")

    command2 = ("java",
                "-Xmx{}g".format(max_mem),
                "-jar",
                "{}".format(config['picard']),
                "BuildBamIndex",
                "INPUT={}".format(output_bam))

    job.fileStore.logToMaster("GATK AddOrReplaceReadGroupsCommand Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    job.fileStore.logToMaster("GATK BuildBamIndex Command: {}\n".format(command2))
    pipeline.run_and_log_command(" ".join(command2), index_log)

    return output_bam


def realign_target_creator(job, config, sample, input_bam, max_mem):
    """Identify targets for realignment"""

    targets = "{}.targets.intervals".format(sample)
    targets_log = "{}.targetcreation.log".format(sample)

    command = ("java",
               "-Xmx{}g".format(max_mem),
               "-jar",
               "{}".format(config['gatk']),
               "-T",
               "IndelRealigner",
               "-R",
               "{}".format(config['reference']),
               "-I",
               "{}".format(input_bam),
               "-o",
               "{}".format(targets),
               "-known",
               "{}".format(config['indel1']),
               "-known",
               "{}".format(config['indel2']),
               "-targetIntervals",
               "{}".format(targets),
               "--read_filter",
               "NotPrimaryAlignment")

    job.fileStore.logToMaster("GATK RealignerTargetCreator Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), targets_log)

    return targets


def realign_indels(job, config, sample, input_bam, targets, max_mem):
    """Create Indel realignment targets and run realignment step"""

    output_bam = "{}.realigned.sorted.bam".format(sample)
    realign_log = "{}.realignindels.log".format(sample)

    command = ("java",
               "-Xmx{}g".format(max_mem),
               "-jar",
               "{}".format(config['gatk']),
               "-T",
               "RealignerTargetCreator",
               "-R",
               "{}".format(config['reference']),
               "-I",
               "{}".format(input_bam),
               "-o",
               "{}".format(targets),
               "-known",
               "{}".format(config['indel1']),
               "-known",
               "{}".format(config['indel2']),
               "-nt",
               "{}".format(multiprocessing.cpu_count()))

    job.fileStore.logToMaster("GATK IndelRealigner Command: {}\n".format(command))
    utilities.touch("{}".format(output_bam))
    pipeline.run_and_log_command(" ".join(command), realign_log)

    return output_bam


def recalibrator(job, config, sample, input_bam, max_mem):
    """Recalibrate and print bases"""

    output_bam = "{}.recalibrated.sorted.bam".format(sample)
    recal_config = "{}.recal".format(sample)
    recal_log = "{}.recalibrate.log".format(sample)
    print_log = "{}.printrecalibrated.log".format(sample)
    cp_log = "{}.copy.log".format(sample)

    # Calculate covariates
    recal_commands = ("java",
                      "-Xmx{}g".format(max_mem),
                      "-jar",
                      "{}".format(config['gatk']),
                      "-T",
                      "BaseRecalibrator",
                      "-R",
                      "{}".format(config['reference']),
                      "-I",
                      "{}".format(input_bam),
                      "-o",
                      "{}".format(recal_config),
                      "--knownSites",
                      "{}".format(config['dbsnp']),
                      "-nct",
                      "{}".format(multiprocessing.cpu_count()))

    # Print recalibrated BAM
    print_reads_command = ("java",
                           "-Xmx{}g".format(max_mem),
                           "-jar",
                           "{}".format(config['gatk']),
                           "-T",
                           "PrintReads",
                           "-R",
                           "{}".format(config['reference']),
                           "-I",
                           "{}".format(input_bam),
                           "-o",
                           "{}".format(output_bam),
                           "-BQSR",
                           "{}".format(recal_config),
                           "-nct",
                           "{}".format(multiprocessing.cpu_count()))

    # Copy index to alternative name
    cp_command = ("cp",
                  "{}.recalibrated.sorted.bai".format(sample),
                  "{}.recalibrated.sorted.bam.bai".format(sample))

    job.fileStore.logToMaster("GATK BaseRecalibrator Command: {}\n".format(recal_commands))
    pipeline.run_and_log_command(" ".join(recal_commands), recal_log)

    job.fileStore.logToMaster("GATK PrintReads Command: {}\n".format(print_reads_command))
    pipeline.run_and_log_command(" ".join(print_reads_command), print_log)

    job.fileStore.logToMaster("GATK Copy Command: {}\n".format(cp_command))
    pipeline.run_and_log_command(" ".join(cp_command), cp_log)

    return output_bam
