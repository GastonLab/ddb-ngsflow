__author__ = 'dgaston'

import sys

from helenus.workflow import pipeline as pipe
from helenus.workflow.elements.utilities import bgzip_and_tabix_vcf_instructions, bgzip_and_tabix_vcf


def run_pindel(project_config, sample_config, tool_config, resource_config):
    """Run pindel on samples"""

    instructions = list()

    # Right now this is inappropriately using parallelization twice. Per sample and also for the total samples
    command_core = ("%s -f %s -c ALL -w %s -E %s -T %s" % (tool_config['pindel']['bin'],
                                                           resource_config['reference_genome'],
                                                           tool_config['pindel']['window'],
                                                           tool_config['pindel']['sensitivity'],
                                                           tool_config['pindel']['num_cores']))

    if project_config['mode'] == "per_sample":
        for sample in sample_config:
            sample['output_prefix'] = "%s_pindel" % sample['name']
            sample['pindel_bam_config'] = "%s.pindel_bam_config.txt" % sample['name']
            logfile = "%s.pindel.log" % sample['name']

            with open(sample['pindel_bam_config'], 'w') as bam_config:
                bam_config.write("%s %s %s\n" % (sample['working_bam'], project_config['insert_size'], sample['name']))

            command = ("%s -o %s -i %s" % (command_core, sample['output_prefix'], sample['pindel_bam_config']))
            sys.stdout.write("Running pindel for sample %s\n" % sample['name'])
            code = pipe.run_and_log_command(command, logfile)
            pipe.check_return_code(code)
    elif project_config['mode'] == "per_cohort":
        sys.stderr.write("ERROR: Not yet implemented\n")
        sys.exit()
    else:
        sys.stderr.write("ERROR: Mode: %s not supported\n" % project_config['mode'])
        sys.exit()

    sys.stdout.write("Finished pindel\n")


def run_pindel2vcf(project_config, sample_config, tool_config, resource_config):
    """Convert pindel output files and convert to VCF, use vcf-convert in stream to convert to vcf version 4.1"""

    instructions = list()

    command_core = ("%s -r %s -R %s -d %s -he 0.01 -G" %
                    (tool_config['pindel2vcf']['bin'], resource_config['reference_genome'],
                     resource_config['genome_build'], resource_config['genome_build']))

    if project_config['mode'] == "per_sample":
        for sample in sample_config:
            sample['pindel_vcf'] = "%s.pindel.vcf" % sample['name']
            if not project_config['ensemble_calling']:
                sample['raw_vcf'] = sample['pindel_vcf']
            logfile = "%s.pindel2vcf.log" % sample['name']

            command = ("%s -P %s -v %s" % (command_core, sample['output_prefix'], sample['pindel_vcf']))

            sample['working_vcf'] = sample['pindel_vcf']
            instructions.append((command, logfile))
    elif project_config['mode'] == "per_cohort":
        sys.stderr.write("ERROR: Mode: %s not supported\n" % project_config['mode'])
        sys.exit()
    else:
        sys.stderr.write("ERROR: Mode: %s not supported\n" % project_config['mode'])
        sys.exit()

    sys.stdout.write("Running pindel2vcf\n")
    pipe.execute_multiprocess(instructions, int(tool_config['pindel2vcf']['num_cores']))
    sys.stdout.write("Finished pindel2vcf\n")