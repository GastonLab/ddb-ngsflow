__author__ = 'dgaston'

import sys

from helenus.workflow import pipeline as pipe
from helenus.workflow.elements.utilities import bgzip_and_tabix_vcf_instructions, bgzip_and_tabix_vcf


def run_varscan_matched():
    """Run VarScan2 in matched tumor/normal mode"""

    raise NotImplementedError()


def run_varscan_unmatched(project_config, sample_config, tool_config, resource_config):
    """Run VarScan2 without a matched normal sample"""

    instructions = list()

    varscan_core = ("java %s -jar %s" % (tool_config['varscan']['max_mem'], tool_config['varscan']['bin']))
    varscan_snp_core = ("%s pileup2snp" % varscan_core)
    varscan_indel_core = ("%s pileup2indel" % varscan_core)
    varscan_cns_core = ("%s pileup2cns" % varscan_core)

    if project_config['mode'] is "per_sample":
        for sample in sample_config:
            logfile = "%s.varscan.log"
            sample['varscan_vcf'] = "%s.varscan.vcf" % sample['name']
            if not project_config['ensemble_calling']:
                sample['raw_vcf'] = sample['varscan_vcf']

            command = ("%s mpileup -f %s %s | %s" % (tool_config['samtools']['bin'],
                                                     resource_config['reference_genome'], sample['working_bam'],
                                                     varscan_core))
            instructions.append((command, logfile))
    elif project_config['mode'] is "per_cohort":
        sys.stderr.write("ERROR: Not yet implemented\n")
        sys.exit()
    else:
        sys.stderr.write("ERROR: Mode: %s not supported\n" % project_config['mode'])

    sys.stdout.write("Running VarDict for samples in project %s\n" % project_config['project_name'])
    pipe.execute_multiprocess(instructions, int(tool_config['varscan']['num_cores']))
    sys.stdout.write("Finished Running VarDict\n")
