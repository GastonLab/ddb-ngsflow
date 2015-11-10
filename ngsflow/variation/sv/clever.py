__author__ = 'dgaston'

import sys

from helenus.workflow import pipeline as pipe
from helenus.workflow.elements.utilities import bgzip_and_tabix_vcf_instructions, bgzip_and_tabix_vcf


def run_clever(project_config, sample_config, tool_config, resource_config):
    """Run CLEVER on samples"""

    instructions = list()
    bgzip_instructions = list()
    tabix_instructions = list()

    for sample in sample_config:
        output_dir = "./%s_CLEVER" % sample['name']
        logfile = "%s.clever.log" % sample['name']
        command = ("%s --use_xa --sorted %s %s %s" %
                   (tool_config['clever']['bin'], sample['working_bam'], resource_config['reference_genome'],
                    output_dir))

        bgzip, tabix = bgzip_and_tabix_vcf_instructions(sample['varscan_vcf'])
        bgzip_instructions.append((bgzip[0], bgzip[1]))
        tabix_instructions.append((tabix[0], tabix[1]))
        instructions.append((command, logfile))

    sys.stdout.write("Running CLEVER for samples in project %s\n" % project_config['project_name'])
    pipe.execute_multiprocess(instructions, int(tool_config['varscan']['num_cores']))
    pipe.execute_multiprocess(bgzip_instructions, int(tool_config['varscan']['num_cores']))
    pipe.execute_multiprocess(tabix_instructions, int(tool_config['varscan']['num_cores']))
    sys.stdout.write("Finished Running CLEVER\n")


def run_laser(project_config, sample_config, tool_config, resource_config):
    """Run LASER on samples"""

    # laser reference.fasta(.gz) reads.1.fastq.gz reads.2.fastq.gz outprefix

    instructions = list()
    bgzip_instructions = list()
    tabix_instructions = list()

    for sample in sample_config:
        output_root = "%s_LASER" % sample['name']
        logfile = "%s.laser.log" % sample['name']
        command = ("%s %s %s %s %s" %
                   (tool_config['laser']['bin'], resource_config['reference_genome'], sample['fastq1'],
                    sample['fastq2'], output_root))

        bgzip, tabix = bgzip_and_tabix_vcf_instructions(sample['varscan_vcf'])
        bgzip_instructions.append((bgzip[0], bgzip[1]))
        tabix_instructions.append((tabix[0], tabix[1]))
        instructions.append((command, logfile))

    sys.stdout.write("Running LASER for samples in project %s\n" % project_config['project_name'])
    pipe.execute_multiprocess(instructions, int(tool_config['varscan']['num_cores']))
    pipe.execute_multiprocess(bgzip_instructions, int(tool_config['varscan']['num_cores']))
    pipe.execute_multiprocess(tabix_instructions, int(tool_config['varscan']['num_cores']))
    sys.stdout.write("Finished Running LASER\n")
