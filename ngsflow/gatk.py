__author__ = 'dgaston'

# This file contains methods for setting up and calling GATK and Picard tool executables in both per_sample and
# per_cohort modes. Haplotype caller, UnifiedGenotyper and MuTect are located in the variation.py file

import sys


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


def run_annotations_and_filters(project_config, sample_config, tool_config, resource_config):
    """GATK Annotate and Variant Filters"""

    instructions1 = list()
    instructions2 = list()

    # Double use of parallelism
    annotation_core = ("java -Xmx%sg -jar %s -T VariantAnnotator -R %s --group StandardAnnotation --dbsnp %s" %
                       (tool_config['gatk']['max_mem'], tool_config['gatk']['bin'],
                        resource_config['reference_genome'], resource_config['dbsnp'],))
    filter_core = ("java -Xmx%sg -jar %s -T VariantFiltration -R %s --filterExpression 'MQ0 > 50' "
                   "--filterName 'HighMQ0' --filterExpression 'DP < 10' --filterName 'LowDepth' "
                   "--filterExpression 'QUAL < 10' --filterName 'LowQual' --filterExpression 'MQ < 10' "
                   "--filterName 'LowMappingQual'" %
                   (tool_config['gatk']['max_mem'], tool_config['gatk']['bin'],
                    resource_config['reference_genome']))

    if project_config['mode'] == "per_sample":
        for sample in sample_config:
            sample['annotated_vcf'] = "%s.annotated.vcf" % sample['name']
            sample['filtered_vcf'] = "%s.filtered.vcf" % sample['name']

            logfile1 = "%s.variantannotation.log" % sample['name']
            logfile2 = "%s.variantfiltration.log" % sample['name']

            command1 = ("%s %s -o %s --variant %s -L %s" % (annotation_core, sample['working_bam'],
                                                            sample['annotated_vcf'], sample['working_vcf'],
                                                            sample['working_vcf']))
            command2 = ("%s -o %s --variant %s" % (filter_core, sample['filtered_vcf'], sample['annotated_vcf']))

            instructions1.append((command1, logfile1))
            instructions2.append((command2, logfile2))
            sample['working_vcf'] = sample['filtered_vcf']
    elif project_config['mode'] == "per_cohort":
        sample_inputs = list()
        for sample in sample_config:
            sample_bam = "-I %s" % sample['working_bam']
            sample_inputs.append(sample_bam)

        sample_bam_string = " ".join(sample_inputs)
        project_config['annotated_vcf'] = "%s.annotated.vcf" % project_config['project_name']
        project_config['filtered_vcf'] = "%s.filtered.vcf" % project_config['project_name']

        logfile1 = "%s.variantannotation.log" % project_config['project_name']
        logfile2 = "%s.variantfiltration.log" % project_config['project_name']

        command1 = ("%s %s -o %s --variant %s -L %s -nt %s" % (annotation_core, sample_bam_string,
                                                               project_config['annotated_vcf'],
                                                               project_config['working_vcf'],
                                                               project_config['working_vcf'],
                                                               tool_config['gatk']['num_cores']))

        command2 = ("%s -o %s --variant %s" % (filter_core, project_config['filtered_vcf'],
                                               project_config['annotated_vcf']))

        instructions1.append((command1, logfile1))
        instructions2.append((command2, logfile2))
        project_config['working_vcf'] = project_config['filtered_vcf']
    else:
        sys.stderr.write("ERROR: Mode: %s not supported\n" % project_config['mode'])
        sys.exit()

    sys.stdout.write("Annotating variants\n")
    pipe.execute_multiprocess(instructions1, int(tool_config['gatk']['num_cores']))

    sys.stdout.write("Applying variant filters\n")
    pipe.execute_multiprocess(instructions2, int(tool_config['gatk']['num_cores']))
    sys.stdout.write("Finished annotating and filtering variants using the GATK\n")


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

    sys.stdout.write("Running MarkDuplicates\n")
    pipe.execute_multiprocess(instructions, int(tool_config['picard']['num_cores']))
    sys.stdout.write("Finished Running MarkDuplicates\n")


def add_or_replace_readgroups(job, config, sample, input_bam):
    """Run AddOrReplaceReadGroups"""

    job.fileStore.logToMaster("Running AddOrReplaceReadGroups in sample: {}".format(sample))

    output_bam = "{}.rg.sorted.bam".format(sample)

    command = ("java",
               "-Xmx{}g".format(config['max_mem']),
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
                "-Xmx{}g".format(config['max_mem']),
                "-jar",
                "{}".format(config['picard']),
                "BuildBamIndex",
                "INPUT={}".format(output_bam))

    job.fileStore.logToMaster("GATK AddOrReplaceReadGroupsCommand Command: {}\n".format(command))
    job.fileStore.logToMaster("GATK BuildBamIndex Command: {}\n".format(command2))

    # p = sub.Popen(command, stdout=sub.PIPE, stderr=err, shell=True)
    # output = p.communicate()
    # code = p.returncode
    # if code:
    #     sys.stdout.write("An error occurred. Please check %s for details\n" % logfile)
    #     sys.stdout.write("%s\n" % output)
    #     sys.stderr.write("An error occurred. Please check %s for details\n" % logfile)

    # p = sub.Popen(command, stdout=sub.PIPE, stderr=err, shell=True)
    # output = p.communicate()
    # code = p.returncode
    # if code:
    #     sys.stdout.write("An error occurred. Please check %s for details\n" % logfile)
    #     sys.stdout.write("%s\n" % output)
    #     sys.stderr.write("An error occurred. Please check %s for details\n" % logfile)

    return output_bam


def realign_indels(job, config, sample, input_bam):
    """Create Indel realignment targets and run realignment step"""

    targets = "{}.targets.intervals".format(sample)
    output_bam = "{}.realigned.sorted.bam".format(sample)

    command = ("java",
               "-Xmx{}g".format(config['max_mem']),
               "-jar",
               "{}".format(config['gatk']),
               "-T",
               "RealignerTargetCreator",
               "-R",
               "{}".format(config['reference']),
               "-known",
               "{}".format(config['indel1']),
               "-known",
               "{}".format(config['indel2']),
               "-I",
               "{}".format(input_bam),
               "-o",
               "{}".format(targets))

    command2 = ("java",
                "-Xmx{}g".format(config['max_mem']),
                "-jar",
                "{}".format(config['gatk']),
                "-T",
                "IndelRealigner",
                "-I",
                "{}".format(input_bam),
                "-o",
                "{}".format(output_bam),
                "-known",
                "{}".format(config['indel1']),
                "-known",
                "{}".format(config['indel2']),
                "-targetIntervals",
                "{}".format(targets),
                "-R",
                "{}".format(config['reference']),
                "--read_filter",
                "NotPrimaryAlignment")

    job.fileStore.logToMaster("GATK RealignerTargetCreator Command: {}\n".format(command))
    job.fileStore.logToMaster("GATK IndelRealigner Command: {}\n".format(command2))

    # p = sub.Popen(command, stdout=sub.PIPE, stderr=err, shell=True)
    # output = p.communicate()
    # code = p.returncode
    # if code:
    #     sys.stdout.write("An error occurred. Please check %s for details\n" % logfile)
    #     sys.stdout.write("%s\n" % output)
    #     sys.stderr.write("An error occurred. Please check %s for details\n" % logfile)

    # p = sub.Popen(command, stdout=sub.PIPE, stderr=err, shell=True)
    # output = p.communicate()
    # code = p.returncode
    # if code:
    #     sys.stdout.write("An error occurred. Please check %s for details\n" % logfile)
    #     sys.stdout.write("%s\n" % output)
    #     sys.stderr.write("An error occurred. Please check %s for details\n" % logfile)




def run_recalibrator(project_config, sample_config, tool_config, resource_config):
    """Recalibrate and print bases"""

    sys.stdout.write("Recalibrating bases for all samples\n")

    instructions1 = list()
    # instructions2 = list()
    instructions3 = list()
    # instructions4 = list()
    instructions5 = list()

    for sample in sample_config:
        sys.stdout.write("Generating commands for multiprocessing steps %s\n" % sample['name'])

        logfile1 = "%s.baserecalibrator.log" % sample['name']
        # logfile2 = "%s.baserecalibrator_second_pass.log" % sample['name']
        logfile3 = "%s.printreads.log" % sample['name']
        # logfile4 = "%s.printplots.log" % sample['name']
        logfile5 = "%s.cpindex.log" % sample['name']

        recal_config = "%s.recal" % sample['name']
        # post_recal= "%s.post.recal" % sample['name']
        # plots = "%s.recalibration_plots.pdf" % sample['name']
        sample['recalibrated_bam'] = "%s.recalibrated.sorted.bam" % sample['name']

        # Calculate covariates
        command1 = ("java -Xmx%sg -jar %s -T BaseRecalibrator -I %s -o %s -R %s --knownSites %s"
                    % (tool_config['gatk']['max_mem'], tool_config['gatk']['bin'], sample['working_bam'],
                       recal_config, resource_config['reference_genome'], resource_config['dbsnp']))

        # Second pass after bqsr
        # command2 = ("java %s -jar %s -T BaseRecalibrator -I %s -o %s -R %s --knownSites %s -BQSR %s"
        #            % (tool_config['gatk']['max_mem'], tool_config['gatk']['bin'], realigned, post_recal,
        # resource_config['reference_genome'], resource_config['dbsnp'], recal_config))

        # Print recalibrated BAM
        command3 = ("java -Xmx%sg -jar %s -T PrintReads -I %s -o %s -R %s -BQSR %s"
                    % (tool_config['gatk']['max_mem'], tool_config['gatk']['bin'], sample['working_bam'],
                       sample['recalibrated_bam'], resource_config['reference_genome'], recal_config))

        # Analysis of Covariates and Plot Printing
        # command4 = ("java %s -jar %s -T AnalyzeCovariates -before %s -after %s -plots %s"
        #            % (tool_config['gatk']['max_mem'], tool_config['gatk']['bin'],recal_config, post_recal,
        #               plots))

        # Copy index to alternative name
        command5 = ("cp %s.recalibrated.sorted.bai %s.recalibrated.sorted.bam.bai" % (sample['name'], sample['name']))

        instructions1.append((command1, logfile1))
        # instructions2.append((command2, logfile2))
        instructions3.append((command3, logfile3))
        # instructions4.append((command4, logfile4))
        instructions5.append((command5, logfile5))

        sample['working_bam'] = sample['recalibrated_bam']

    sys.stdout.write("Running multiprocessing of BaseRecalibrator\n")
    pipe.execute_multiprocess(instructions1, int(tool_config['gatk']['num_cores']))

    #
    # sys.stdout.write("Running multiprocessing of Post Calibration BaseRecalibrator\n")
    # pipe.execute_multiprocess(instructions2, int(tool_config['gatk']['num_cores']))

    sys.stdout.write("Running multiprocessing of PrintReads\n")
    pipe.execute_multiprocess(instructions3, int(tool_config['gatk']['num_cores']))

    # sys.stdout.write("Running Analysis and Plotting of Covariates\n")
    # pipe.execute_multiprocess(instructions4, int(tool_config['gatk']['num_cores']))

    sys.stdout.write("Copying index files to alternative name for special tools\n")
    pipe.execute_multiprocess(instructions5, tool_config['gatk']['num_cores'])

    sys.stdout.write("Finished recalibrating bases\n")


