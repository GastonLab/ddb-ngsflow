__author__ = 'dgaston'

import sys

from helenus.workflow import pipeline as pipe
from helenus.workflow.elements.utilities import bgzip_and_tabix_vcf_instructions, bgzip_and_tabix_vcf


def merge_variant_calls(project_config, sample_config, tool_config, resource_config):
    """Use vcf-isec to merge vcfs together and create an ensemble call set report"""

    instructions = list()
    bgzip_instructions = list()
    tabix_instructions = list()

    isec_command_core = ("%s -n +1" % (tool_config['vcftools_isec']['bin']))

    if project_config['mode'] == "per_sample":
        for sample in sample_config:
            if 'variant_callers' in tool_config.keys():
                vcf_files = list()
                for caller in tool_config['variant_callers']:
                    bgzip, tabix = bgzip_and_tabix_vcf_instructions(sample[caller])
                    bgzip_instructions.append((bgzip[0], bgzip[1]))
                    tabix_instructions.append((tabix[0], tabix[1]))
                    vcf_files.append("%s.gz" % sample[caller])
                vcf_files_string = " ".join(vcf_files)

                logfile = "%s.mergevariants.log" % sample['name']
                sample['merged_vcf'] = "%s.merged.vcf" % sample['name']
                sample['raw_vcf'] = sample['merged_vcf']
                sample['working_vcf'] = sample['merged_vcf']

                isec_command = ("%s %s > %s" % (isec_command_core, vcf_files_string, sample['merged_vcf']))

                instructions.append((isec_command, logfile))
            if 'sv_callers' in tool_config.keys():
                vcf_files = list()
                for caller in tool_config['sv_callers']:
                    bgzip, tabix = bgzip_and_tabix_vcf_instructions(sample[caller])
                    bgzip_instructions.append((bgzip[0], bgzip[1]))
                    tabix_instructions.append((tabix[0], tabix[1]))
                    vcf_files.append("%s.gz" % sample[caller])
                vcf_files_string = " ".join(vcf_files)

                logfile = "%s.mergesvvariants.log" % sample['name']
                sample['merged_sv_vcf'] = "%s.merged-sv.vcf" % sample['name']
                sample['raw_sv_vcf'] = sample['merged_sv_vcf']
                sample['working_sv_vcf'] = sample['merged_sv_vcf']

                isec_command = ("%s %s > %s" % (isec_command_core, vcf_files_string, sample['merged_sv_vcf']))

                instructions.append((isec_command, logfile))
    elif project_config['mode'] == "per_cohort":
        if 'variant_callers' in tool_config.keys():
            vcf_files = list()
            for caller in tool_config['variant_callers']:
                vcf_files.append(project_config[caller])
            vcf_files_string = " ".join(vcf_files)

            logfile = "%s.mergevariants.log" % project_config['project_name']
            project_config['merged_vcf'] = "%s.mergedvariants.vcf" % project_config['project_name']
            project_config['raw_vcf'] = project_config['merged_vcf']
            project_config['working_vcf'] = project_config['merged_vcf']

            isec_command = ("%s %s > %s" % (isec_command_core, vcf_files_string, project_config['merged_vcf']))
            instructions.append((isec_command, logfile))
        if 'sv_callers' in tool_config.keys():
                vcf_files = list()
                for caller in tool_config['sv_callers']:
                    bgzip, tabix = bgzip_and_tabix_vcf_instructions(project_config[caller])
                    bgzip_instructions.append((bgzip[0], bgzip[1]))
                    tabix_instructions.append((tabix[0], tabix[1]))
                    vcf_files.append("%s.gz" % project_config[caller])
                vcf_files_string = " ".join(vcf_files)

                logfile = "%s.mergesvvariants.log" % project_config['project_name']
                project_config['merged_sv_vcf'] = "%s.merged-sv.vcf" % project_config['project_name']
                project_config['raw_sv_vcf'] = project_config['merged_sv_vcf']
                project_config['working_sv_vcf'] = project_config['merged_sv_vcf']

                isec_command = ("%s %s > %s" % (isec_command_core, vcf_files_string, project_config['merged_sv_vcf']))

                instructions.append((isec_command, logfile))
    else:
        sys.stderr.write("ERROR: Mode: %s not supported\n" % project_config['mode'])
        sys.exit()

    sys.stdout.write("Running BGZip and Tabix on all input VCF files prior to merging\n")
    pipe.execute_multiprocess(bgzip_instructions, int(tool_config['vcftools_isec']['num_cores']))
    pipe.execute_multiprocess(tabix_instructions, int(tool_config['vcftools_isec']['num_cores']))

    sys.stdout.write("Running vcftools isec\n")
    pipe.execute_multiprocess(instructions, int(tool_config['vcftools_isec']['num_cores']))
    sys.stdout.write("Finished vcftools isec\n")


def evaluate_double_strand_calls(project_config, sample_config, tool_config, resource_config):
    """Use pool A and pool B data from an Illumina strand-specific sequencing experiment to call high confidences
    and likely false-positive variants"""

    new_samples_config = list()
    instructions = list()
    bgzip_instructions = list()
    tabix_instructions = list()

    isec_command_core = ("%s -n +1 -f" % (tool_config['vcftools_isec']['bin']))

    if project_config['mode'] == "per_sample":
        for sample in project_config['stranded_sample_data']:
            new_sample_dict = dict()
            pool_dicts = list()
            vcf_files = list()

            logfile = "%s.mergepools.log" % sample['name']
            pool1_name = sample['pools'][0]
            pool2_name = sample['pools'][1]

            for sample_dict in sample_config:
                if sample_dict['name'] == pool1_name or sample_dict['name'] == pool2_name:
                    pool_dicts.append(sample_dict)
                    bgzip, tabix = bgzip_and_tabix_vcf_instructions(sample_dict['working_vcf'])
                    bgzip_instructions.append((bgzip[0], bgzip[1]))
                    tabix_instructions.append((tabix[0], tabix[1]))
                    vcf_files.append("%s.gz" % sample_dict['working_vcf'])

            new_sample_dict['name'] = sample['name']
            new_sample_dict['double_strand_merged_vcf'] = "%s.ds.merged.vcf" % sample['name']
            new_sample_dict['raw_vcf'] = new_sample_dict['double_strand_merged_vcf']
            new_sample_dict['working_vcf'] = new_sample_dict['double_strand_merged_vcf']
            new_samples_config.append(new_sample_dict)

            vcf_files_string = " ".join(vcf_files)
            isec_command = ("%s %s > %s" % (isec_command_core, vcf_files_string,
                                            new_sample_dict['double_strand_merged_vcf']))

            instructions.append((isec_command, logfile))

            if 'sv_callers' in tool_config.keys():
                new_sample_dict = dict()
                pool_dicts = list()
                vcf_files = list()

                logfile = "%s.mergesvpools.log" % sample['name']
                pool1_name = sample['pools'][0]
                pool2_name = sample['pools'][1]

                for sample_dict in sample_config:
                    if sample_dict['name'] == pool1_name or sample_dict['name'] == pool2_name:
                        pool_dicts.append(sample_dict)
                        bgzip, tabix = bgzip_and_tabix_vcf_instructions(sample_dict['working_sv_vcf'])
                        bgzip_instructions.append((bgzip[0], bgzip[1]))
                        tabix_instructions.append((tabix[0], tabix[1]))
                        vcf_files.append("%s.gz" % sample_dict['working_sv_vcf'])

                new_sample_dict['name'] = sample['name']
                new_sample_dict['double_strand_merged_sv_vcf'] = "%s.ds.merged.sv.vcf" % sample['name']
                new_sample_dict['raw_sv_vcf'] = new_sample_dict['double_strand_merged_sv_vcf']
                new_sample_dict['working_sv_vcf'] = new_sample_dict['double_strand_merged_sv_vcf']
                new_samples_config.append(new_sample_dict)

                vcf_files_string = " ".join(vcf_files)
                isec_command = ("%s %s > %s" % (isec_command_core, vcf_files_string,
                                                new_sample_dict['double_strand_merged_sv_vcf']))

                instructions.append((isec_command, logfile))

    elif project_config['mode'] == "per_cohort":
        sys.stdout.write("ERROR: Double-Stranded variant call merging not supported for per cohort variants\n")
        sys.exit()
    else:
        sys.stderr.write("ERROR: Mode: %s not supported\n" % project_config['mode'])
        sys.exit()

    sys.stdout.write("Running BGZip and Tabix on all input VCF files prior to merging\n")
    pipe.execute_multiprocess(bgzip_instructions, int(tool_config['vcftools_isec']['num_cores']))
    pipe.execute_multiprocess(tabix_instructions, int(tool_config['vcftools_isec']['num_cores']))

    sys.stdout.write("Running vcftools isec for double stranded pools\n")
    pipe.execute_multiprocess(instructions, int(tool_config['vcftools_isec']['num_cores']))
    sys.stdout.write("Finished vcftools isec\n")

    # Set new configuration parameters for new merged samples
    del sample_config[:]
    for sample in new_samples_config:
        sample_config.append(sample)
