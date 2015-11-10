__author__ = 'dgaston'

import sys

from helenus.workflow import pipeline as pipe
from helenus.workflow.elements.utilities import bgzip_and_tabix_vcf_instructions, bgzip_and_tabix_vcf


def run_platypus_single(project_config, sample_config, tool_config, resource_config):
    """Run platypus on a single sample"""

    instructions = list()

    sys.stdout.write("Running Platypus for samples in project %s\n" % project_config['project_name'])

    for sample in sample_config:
        sample['platypus_vcf'] = "%s.platypus.vcf" % sample['name']

        if not project_config['ensemble_calling']:
            sample['raw_vcf'] = sample['platypus_vcf']
            sample['working_vcf'] = sample['platypus_vcf']

        logfile = "%s.platypus_single.log" % sample['name']
        internal_log = "%s.platypus_internal.log" % sample['name']

        command = ("%s callVariants --refFile=%s --regions=%s --assemble=1 --assembleBrokenPairs=1 "
                   "--filterDuplicates=0 --nCPU=%s --logFileName=%s --bamFiles=%s --output=%s" %
                   (tool_config['platypus']['bin'], resource_config['reference_genome'],
                    tool_config['platypus']['regions'], tool_config['platypus']['num_cores'], internal_log,
                    sample['working_bam'], sample['platypus_vcf']))

        sample['working_vcf'] = sample['platypus_vcf']

        code = pipe.run_and_log_command(command, logfile)
        pipe.check_return_codes(code)
        instructions.append((command, logfile))

    sys.stdout.write("Executing Platypus\n")
    pipe.execute_multiprocess(instructions, int(tool_config['platypus']['num_cores']))
    sys.stdout.write("Finished Running Platypus\n")
