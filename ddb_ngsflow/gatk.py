"""
.. module:: gatk
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling GATK utilities.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>

"""

import pipeline


def diagnosetargets(job, config, name, samples, input_bam):
    """Run GATK's DiagnoseTargets against the supplied region
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param samples: samples dictionary.
    :type samples: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The DiagnoseTargets output vcf file name.
    """

    diagnose_targets_vcf = "{}.diagnosetargets.vcf".format(name)
    missing_intervals = "{}.missing.intervals".format(name)
    logfile = "{}.diagnose_targets.log".format(name)

    command = ("{}".format(config['gatk']['bin']),
               "-T",
               "DiagnoseTargets",
               "-R",
               "{}".format(config['reference']),
               "-L",
               "{}".format(samples[name]['regions']),
               "--coverage_status_threshold",
               "{}".format(config['coverage_loci_threshold']),
               "--bad_mate_status_threshold",
               "{}".format(config['bad_mate_threshold']),
               "--minimum_coverage",
               "{}".format(config['coverage_threshold']),
               "--quality_status_threshold",
               "{}".format(config['quality_loci_threshold']),
               "-I",
               "{}".format(input_bam),
               "-o",
               "{}".format(diagnose_targets_vcf),
               "--missing_intervals",
               "{}".format(missing_intervals))

    job.fileStore.logToMaster("GATK DiagnoseTargets Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return diagnose_targets_vcf


def diagnose_pooled_targets(job, config, name, regions, samples, input_bam1, input_bam2):
    """Run GATK's DiagnoseTargets against the supplied region
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param regions: regions dictionary key name and tag.
    :type regions: str.
    :param samples: samples dictionary.
    :type samples: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The DiagnoseTargets output vcf file name.
    """

    diagnose_targets_vcf = "{}_{}.diagnosetargets.vcf".format(name, regions)
    missing_intervals = "{}_{}.missing.intervals".format(name, regions)
    logfile = "{}.{}.diagnose_targets.log".format(name, regions)

    command = ("{}".format(config['gatk']['bin']),
               "-T",
               "DiagnoseTargets",
               "-R",
               "{}".format(config['reference']),
               "-L",
               "{}".format(samples[name][regions]),
               "--coverage_status_threshold",
               "{}".format(config['coverage_loci_threshold']),
               "--bad_mate_status_threshold",
               "{}".format(config['bad_mate_threshold']),
               "--minimum_coverage",
               "{}".format(config['coverage_threshold']),
               "--quality_status_threshold",
               "{}".format(config['quality_loci_threshold']),
               "-I",
               "{}".format(input_bam1),
               "-I",
               "{}".format(input_bam2),
               "-o",
               "{}".format(diagnose_targets_vcf),
               "--missing_intervals",
               "{}".format(missing_intervals))

    job.fileStore.logToMaster("GATK DiagnoseTargets Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return diagnose_targets_vcf


def annotate_vcf(job, config, name, input_vcf, input_bam):
    """Run GATK's VariantAnnotation on the specified VCF
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_vcf: The input_vcf file name to process.
    :type input_vcf: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    output_vcf = "{}.annotated.vcf".format(name)
    annotation_logfile = "{}.variantannotation.log".format(name)

    annotation_command = ("{}".format(config['gatk']['bin']),
                          "-T",
                          "VariantAnnotator",
                          "-R",
                          "{}".format(config['reference']),
                          "-nt",
                          "{}".format(config['gatk']['num_cores']),
                          "--group",
                          "StandardAnnotation",
                          "--dbsnp",
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


def filter_variants(job, config, name, input_vcf):
    """Run GATK's VariantFilter on the specified VCF
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_vcf: The input_vcf file name to process.
    :type input_vcf: str.
    :returns:  str -- The output vcf file name.
    """

    output_vcf = "{}.filtered.vcf".format(name)
    filter_log = "{}.variantfiltration.log".format(name)

    filter_command = ("{}".format(config['gatk']['bin']),
                      "-T",
                      "VariantFiltration",
                      "-R",
                      "{}".format(config['reference']),
                      "--filterExpression",
                      "'MQ0 > {}'".format(config['mq0_threshold']),
                      "--filterName",
                      "'HighMQ0'",
                      "--filterExpression",
                      "'DP < {}'".format(config['coverage_threshold']),
                      "--filterName",
                      "'LowDepth'",
                      "--filterExpression",
                      "'QUAL < {}'".format(config['var_qual_threshold']),
                      "--filterName",
                      "'LowQual'",
                      "--filterExpression",
                      "'MQ < {}'".format(config['map_qual_threshold']),
                      "--filterName",
                      "'LowMappingQual'",
                      "--variant",
                      "{}".format(input_vcf),
                      "-o",
                      "{}".format(output_vcf))

    job.fileStore.logToMaster("GATK VariantFiltration Command: {}\n".format(filter_command))
    pipeline.run_and_log_command(" ".join(filter_command), filter_log)

    return output_vcf


def mark_duplicates(job, config, name, input_bam):
    """Run Picard MarkDuplicates
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output bam file name.
    """

    job.fileStore.logToMaster("Running MarkDuplicates for sample: {}".format(name))

    metrics_file = "{}.dedup.metrics".format(name)
    output_bam = "{}.dedup.sorted.bam".format(name)
    logfile = "{}.markduplicates.log".format(name)

    command = ("{}".format(config['picard']['bin']),
               "MarkDuplicates",
               "CREATE_INDEX=true",
               "METRICS_FILE={}".format(metrics_file),
               "INPUT={}".format(input_bam),
               "OUTPUT={}".format(output_bam))

    job.fileStore.logToMaster("Picard MarkDuplicates Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return output_bam


def add_or_replace_readgroups(job, config, name, input_bam):
    """Run Picard's AddOrReplaceReadGroups on the specified BAM
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output bam file name.
    """

    job.fileStore.logToMaster("Running AddOrReplaceReadGroups in sample: {}".format(name))

    output_bam = "{}.rg.sorted.bam".format(name)
    logfile = "{}.addreadgroups.log".format(name)
    index_log = "{}.buildindex.log".format(name)

    command = ("{}".format(config['picard']['bin']),
               "AddOrReplaceReadGroups",
               "INPUT={}".format(input_bam),
               "OUTPUT={}".format(output_bam),
               "RGID={}".format(name),
               "RGSM={}".format(name),
               "RGLB={}".format(name),
               "RGPL=illumina",
               "RGPU=miseq")

    command2 = ("{}".format(config['picard']['bin']),
                "BuildBamIndex",
                "INPUT={}".format(output_bam))

    job.fileStore.logToMaster("GATK AddOrReplaceReadGroupsCommand Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    job.fileStore.logToMaster("GATK BuildBamIndex Command: {}\n".format(command2))
    pipeline.run_and_log_command(" ".join(command2), index_log)

    return output_bam


def realign_target_creator(job, config, name, input_bam):
    """Run GATK TargetCreator on the specified BAM to identify targets for realignment
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The file name of the targets file.
    """

    targets = "{}.targets.intervals".format(name)
    targets_log = "{}.targetcreation.log".format(name)

    command = ("{}".format(config['gatk']['bin']),
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
               "{}".format(config['gatk']['num_cores'])
               )

    job.fileStore.logToMaster("GATK RealignerTargetCreator Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), targets_log)

    return targets


def realign_indels(job, config, name, input_bam, targets):
    """Run GATK Indel Realignment on the specified BAM
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :param targets: The file name of targets to realign.
    :type targets: str.
    :returns:  str -- The output bam file name.
    """

    output_bam = "{}.realigned.sorted.bam".format(name)
    realign_log = "{}.realignindels.log".format(name)

    command = ("{}".format(config['gatk']['bin']),
               "-T",
               "IndelRealigner",
               "-R",
               "{}".format(config['reference']),
               "-I",
               "{}".format(input_bam),
               "-known",
               "{}".format(config['indel1']),
               "-known",
               "{}".format(config['indel2']),
               "-targetIntervals",
               "{}".format(targets),
               "--read_filter",
               "NotPrimaryAlignment",
               "-o",
               "{}".format(output_bam))

    job.fileStore.logToMaster("GATK IndelRealigner Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), realign_log)

    return output_bam


def recalibrator(job, config, name, input_bam):
    """Run GATK Recalibrator on the specified BAM

    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output bam file name.

    """

    output_bam = "{}.recalibrated.sorted.bam".format(name)
    recal_config = "{}.recal".format(name)
    recal_log = "{}.recalibrate.log".format(name)
    print_log = "{}.printrecalibrated.log".format(name)
    cp_log = "{}.copy.log".format(name)

    # Calculate covariates
    recal_commands = ("{}".format(config['gatk']['bin']),
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
                      "{}".format(config['gatk']['num_cores']))

    # Print recalibrated BAM
    print_reads_command = ("{}".format(config['gatk']['bin']),
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
                           "{}".format(config['gatk']['num_cores']))

    # Copy index to alternative name
    cp_command = ("cp",
                  "{}.recalibrated.sorted.bai".format(name),
                  "{}.recalibrated.sorted.bam.bai".format(name))

    job.fileStore.logToMaster("GATK BaseRecalibrator Command: {}\n".format(recal_commands))
    pipeline.run_and_log_command(" ".join(recal_commands), recal_log)

    job.fileStore.logToMaster("GATK PrintReads Command: {}\n".format(print_reads_command))
    pipeline.run_and_log_command(" ".join(print_reads_command), print_log)

    job.fileStore.logToMaster("GATK Copy Command: {}\n".format(cp_command))
    pipeline.run_and_log_command(" ".join(cp_command), cp_log)

    return output_bam
