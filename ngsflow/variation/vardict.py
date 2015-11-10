__author__ = 'dgaston'

import sys

from helenus.workflow import pipeline as pipe
from helenus.workflow.elements.utilities import bgzip_and_tabix_vcf_instructions, bgzip_and_tabix_vcf


def run_vardict_matched():
    """Run VarDict in matched tumor/normal mode"""

    raise NotImplementedError()


def run_vardict_single(project_config, sample_config, tool_config, resource_config):
    """Run VarDict without a matched normal sample"""

    instructions = list()

    vardict_core = ("%s -G %s -z -c 1 -S 2 -E 3 -g 4 -f %s" % (tool_config['vardict']['bin'],
                                                               resource_config['reference_genome'],
                                                               tool_config['min_alt_af']))

    vardict2vcf_core = ("%s -E -f %s" % (tool_config['vardict2vcf']['bin'], tool_config['min_alt_af']))

    if project_config['mode'] == "per_sample":
        for sample in sample_config:
            logfile = "%s.vardict.log" % sample['name']
            sample['vardict_vcf'] = "%s.vardict.vcf" % sample['name']
            if not project_config['ensemble_calling']:
                sample['raw_vcf'] = sample['vardict_vcf']

            command = ("%s -N %s -b %s %s | %s | %s -N %s > %s" %
                       (vardict_core, sample['name'], sample['recalibrated_bam'], resource_config['regions'],
                        tool_config['vardict_strandbias']['bin'], vardict2vcf_core, sample['name'],
                        sample['vardict_vcf']))

            instructions.append((command, logfile))
    elif project_config['mode'] == "per_cohort":
        sys.stderr.write("ERROR: Not yet implemented\n")
        sys.exit()
    else:
        sys.stderr.write("ERROR: Mode: %s not supported\n" % project_config['mode'])

    sys.stdout.write("Running VarDict for samples in project %s\n" % project_config['project_name'])
    pipe.execute_multiprocess(instructions, int(tool_config['vardict']['num_cores']))
    sys.stdout.write("Finished Running VarDict\n")