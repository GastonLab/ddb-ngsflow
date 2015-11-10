__author__ = 'dgaston'

import sys

from helenus.workflow import pipeline as pipe
from helenus.workflow.elements.utilities import bgzip_and_tabix_vcf_instructions, bgzip_and_tabix_vcf


def run_delly2_single(project_config, sample_config, tool_config, resource_config):
    """Run delly2 for structural variant detection. As delly2 is parallelized on the level of samples,
    we use a single-threaded version"""

    instructions = list()
    bgzip_instructions = list()
    tabix_instructions = list()

    delly_command_core = ("%s -x %s -g %s" % (tool_config['delly']['bin'], resource_config['delly']['exclude'],
                                              resource_config['reference_genome']))

    if project_config['mode'] is "per_sample":
        for sample in sample_config:
            sample['delly_vcfs'] = list()
            for mut_type in ["DEL", "DUP", "TRA", "INV"]:
                output_vcf = "%s.%s.vcf" % (sample['name'], mut_type)
                logfile = "%s.%s.log" % (sample['name'], mut_type)
                sample['delly_vcfs'].append(output_vcf)
                command = ("%s -t %s -o %s %s" % (delly_command_core, mut_type, output_vcf, sample['working_bam']))

                bgzip, tabix = bgzip_and_tabix_vcf_instructions(output_vcf)
                bgzip_instructions.append((bgzip[0], bgzip[1]))
                tabix_instructions.append((tabix[0], tabix[1]))
                instructions.append((command, logfile))
    elif project_config['mode'] is "per_cohort":
        sample_bams = list()
        for sample in sample_config:
            sample_bams.append(sample['working_bam'])
        sample_bams_list = " ".join(sample_bams)
        project_config['delly_vcfs'] = list()
        for mut_type in ["DEL", "DUP", "TRA", "INV"]:
            output_vcf = "%s.%s.vcf" % (project_config['project_name'], mut_type)
            logfile = "%s.%s.log" % (project_config['project_name'], mut_type)
            project_config['delly_vcfs'].append(output_vcf)
            command = ("%s -t %s -o %s %s" % (delly_command_core, mut_type, output_vcf, sample_bams_list))

            bgzip, tabix = bgzip_and_tabix_vcf_instructions(output_vcf)
            bgzip_instructions.append((bgzip[0], bgzip[1]))
            tabix_instructions.append((tabix[0], tabix[1]))
            instructions.append((command, logfile))
    else:
        sys.stderr.write("ERROR: Mode: %s not supported\n" % project_config['mode'])
        sys.exit()

    sys.stdout.write("Running Delly\n")
    pipe.execute_multiprocess(instructions, int(tool_config['delly']['num_cores']))
    pipe.execute_multiprocess(bgzip_instructions, int(tool_config['delly']['num_cores']))
    pipe.execute_multiprocess(tabix_instructions, int(tool_config['delly']['num_cores']))
    sys.stdout.write("Finished Delly\n")

    sys.stdout.write("Merging Delly output\n")

    sys.stdout.write("Finished merging Delly output\n")
