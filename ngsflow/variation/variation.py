__author__ = 'dgaston'

import sys
import time
import multiprocessing
import subprocess as sub

from ngsflow.utils import utilities
from .. import pipeline as pipe


def merge_variant_calls(job, config, sample, vcf_files):
    """Use vcf-isec to merge vcfs together and create an ensemble call set report"""

    for vcf in vcf_files:
        utilities.bgzip_and_tabix_vcf(vcf)
        vcf_files.append("{}.gz".format(vcf))
    vcf_files_string = " ".join(vcf_files)

    merged_vcf = "{}.merged.vcf".format(sample)
    logfile = "{}.merged.log".format(sample)

    isec_command = ("vcf-isec",
                    "-n",
                    "+1",
                    "{}".format(vcf_files_string),
                    ">",
                    "{}".format(merged_vcf))

    job.fileStore.logToMaster("Vcftools intersect Command: {}\n".format(isec_command))
    # pipeline.run_and_log_command(" ".join(isec_command), logfile)
    utilities.touch(merged_vcf)
    time.sleep(2)

    return merged_vcf


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
