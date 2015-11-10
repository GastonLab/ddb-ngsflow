__author__ = 'dgaston'

import sys

from helenus.workflow import pipeline as pipe
from helenus.workflow.elements.utilities import bgzip_and_tabix_vcf_instructions, bgzip_and_tabix_vcf


# Indelminer doesn't have built-in multithreading or parallelization?
def run_indelminer_single(project_config, sample_config, tool_config, resource_config):
    """Run indelminer for calling indels in sample data"""

    instructions = list()

    indelminer_command_core = ("%s %s" % (tool_config['indelminer']['bin'], resource_config['reference_genome']))

    for sample in sample_config:
        sample['indelminer_vcf'] = "%s.indelminer.vcf" % sample['name']
        if not project_config['ensemble_calling']:
            sample['raw_vcf'] = sample['indelminer_vcf']
        logfile = "%s.indelminer.log" % sample['name']

        command = ("%s sample=%s > %s" % (indelminer_command_core, sample['working_bam'], sample['indelminer_vcf']))
        instructions.append((command, logfile))

    sys.stdout.write("Running Indelminer\n")
    pipe.execute_multiprocess(instructions, int(tool_config['indelminer']['num_cores']))
    sys.stdout.write("Finished Indelminer\n")
