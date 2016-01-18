"""
.. module:: gatk
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling GATK utilities.

.. moduleauthor:: Daniel Gaston <daniel.gaston@gmail.com>

"""

from ngsflow import pipeline


def diagnosetargets(job, config, sample, samples, input_bam):
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

    diagnose_targets_vcf = "{}.diagnosetargets.vcf".format(sample)
    missing_intervals = "{}.missing.intervals".format(sample)
    logfile = "{}.diagnose_targets.log".format(sample)

    command = ("java",
               "-Xmx{}g".format(config['gatk']['max_mem']),
               "-jar",
               "{}".format(config['gatk']['bin']),
               "-T",
               "DiagnoseTargets",
               "-R",
               "{}".format(config['reference']),
               "-L",
               "{}".format(samples[sample]['regions']),
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


def annotate_vcf(job, config, sample, input_vcf, input_bam):
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

    output_vcf = "{}.annotated.vcf".format(sample)
    annotation_logfile = "{}.variantannotation.log".format(sample)

    annotation_command = ("java",
                          "-Xmx{}g".format(config['gatk']['max_mem']),
                          "-jar",
                          "{}".format(config['gatk']['bin']),
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


def filter_variants(job, config, sample, input_vcf):
    """Run GATK's VariantFilter on the specified VCF
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_vcf: The input_vcf file name to process.
    :type input_vcf: str.
    :returns:  str -- The output vcf file name.
    """

    output_vcf = "{}.filtered.vcf".format(sample)
    filter_log = "{}.variantfiltration.log".format(sample)

    filter_command = ("java",
                      "-Xmx{}g".format(config['gatk']['max_mem']),
                      "-jar",
                      "{}".format(config['gatk']['bin']),
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


def run_mark_duplicates(job, config, sample, input_bam):
    """Run Picard MarkDuplicates
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output bam file name.
    """

    job.fileStore.logToMaster("Running MarkDuplicates for sample: {}".format(sample))

    metrics_file = "{}.dedup.metrics".format(sample)
    output_bam = "{}.dedup.sorted.bam".format(sample)
    logfile = "{}.markduplicates.log".format(sample)

    command = ("java",
               "-Xmx{}g".format(config['gatk']['max_mem']),
               "-jar",
               "{}".format(config['picard']['bin']),
               "MarkDuplicates",
               "CREATE_INDEX=true",
               "METRICS_FILE={}".format(metrics_file),
               "INPUT={}".format(input_bam),
               "OUTPUT={}".format(output_bam))

    job.fileStore.logToMaster("GATK BuildBamIndex Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return output_bam


def add_or_replace_readgroups(job, config, sample, input_bam):
    """Run Picard's AddOrReplaceReadGroups on the specified BAM
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output bam file name.
    """

    job.fileStore.logToMaster("Running AddOrReplaceReadGroups in sample: {}".format(sample))

    output_bam = "{}.rg.sorted.bam".format(sample)
    logfile = "{}.addreadgroups.log".format(sample)
    index_log = "{}.buildindex.log".format(sample)

    command = ("java",
               "-Xmx{}g".format(config['gatk']['max_mem']),
               "-jar",
               "{}".format(config['picard']['bin']),
               "AddOrReplaceReadGroups",
               "INPUT={}".format(input_bam),
               "OUTPUT={}".format(output_bam),
               "RGID={}".format(sample),
               "RGSM={}".format(sample),
               "RGLB={}".format(sample),
               "RGPL=illumina",
               "RGPU=miseq")

    command2 = ("java",
                "-Xmx{}g".format(config['gatk']['max_mem']),
                "-jar",
                "{}".format(config['picard']['bin']),
                "BuildBamIndex",
                "INPUT={}".format(output_bam))

    job.fileStore.logToMaster("GATK AddOrReplaceReadGroupsCommand Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    job.fileStore.logToMaster("GATK BuildBamIndex Command: {}\n".format(command2))
    pipeline.run_and_log_command(" ".join(command2), index_log)

    return output_bam


def realign_target_creator(job, config, sample, input_bam):
    """Run GATK TargetCreator on the specified BAM to identify targets for realignment
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The file name of the targets file.
    """

    targets = "{}.targets.intervals".format(sample)
    targets_log = "{}.targetcreation.log".format(sample)

    command = ("java",
               "-Xmx{}g".format(config['gatk']['max_mem']),
               "-jar",
               "{}".format(config['gatk']['bin']),
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


def realign_indels(job, config, sample, input_bam, targets):
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

    output_bam = "{}.realigned.sorted.bam".format(sample)
    realign_log = "{}.realignindels.log".format(sample)

    command = ("java",
               "-Xmx{}g".format(config['gatk']['max_mem']),
               "-jar",
               "{}".format(config['gatk']['bin']),
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


def recalibrator(job, config, sample, input_bam):
    """Run GATK Recalibrator on the specified BAM

    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output bam file name.

    """

    output_bam = "{}.recalibrated.sorted.bam".format(sample)
    recal_config = "{}.recal".format(sample)
    recal_log = "{}.recalibrate.log".format(sample)
    print_log = "{}.printrecalibrated.log".format(sample)
    cp_log = "{}.copy.log".format(sample)

    # Calculate covariates
    recal_commands = ("java",
                      "-Xmx{}g".format(config['gatk']['max_mem']),
                      "-jar",
                      "{}".format(config['gatk']['bin']),
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
    print_reads_command = ("java",
                           "-Xmx{}g".format(config['gatk']['max_mem']),
                           "-jar",
                           "{}".format(config['gatk']['bin']),
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
                  "{}.recalibrated.sorted.bai".format(sample),
                  "{}.recalibrated.sorted.bam.bai".format(sample))

    job.fileStore.logToMaster("GATK BaseRecalibrator Command: {}\n".format(recal_commands))
    pipeline.run_and_log_command(" ".join(recal_commands), recal_log)

    job.fileStore.logToMaster("GATK PrintReads Command: {}\n".format(print_reads_command))
    pipeline.run_and_log_command(" ".join(print_reads_command), print_log)

    job.fileStore.logToMaster("GATK Copy Command: {}\n".format(cp_command))
    pipeline.run_and_log_command(" ".join(cp_command), cp_log)

    return output_bam
