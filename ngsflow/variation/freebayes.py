__author__ = 'dgaston'

import sys

from helenus.workflow import pipeline as pipe
from helenus.workflow.elements.utilities import bgzip_and_tabix_vcf_instructions, bgzip_and_tabix_vcf


def run_freebayes_single(project_config, sample_config, tool_config, resource_config):
    """Run FreeBayes without a matched normal sample"""

    instructions = list()

    command_core = ("%s --fasta-reference %s -t %s --min-alternate-fraction %s --pooled-discrete "
                    "--pooled-continuous --genotype-qualities --report-genotype-likelihood-max "
                    "--allele-balance-priors-off --use-duplicate-reads --min-repeat-entropy 1" %
                    (tool_config['freebayes']['bin'], resource_config['reference_genome'],
                     resource_config['merged_regions'], tool_config['min_alt_af']))

    if project_config['mode'] == "per_sample":
        for sample in sample_config:
            sample['freebayes_vcf'] = "%s.freebayes.vcf" % sample['name']

            if not project_config['ensemble_calling']:
                sample['raw_vcf'] = sample['freebayes_vcf']

            logfile = "%s.freebayes.log" % sample['name']

            command = ("%s -v %s %s" % (command_core, sample['freebayes_vcf'], sample['working_bam']))

            sample['working_vcf'] = sample['freebayes_vcf']
            instructions.append((command, logfile))
    elif project_config == "per_cohort":
        bam_file_list = "freebayes_bam_list.txt"
        with open(bam_file_list, 'w') as outfile:
            for sample in sample_config:
                outfile.write("%s\n" % sample['working_bam'])

        project_config['freebayes_vcf'] = "%s.freebayes.vcf" % project_config['project_name']

        if not project_config['ensemble_calling']:
            project_config['raw_vcf'] = project_config['freebayes_vcf']

        logfile = "%s.freebayes.log" % project_config['project_name']

        command = ("%s  -v %s %s" % (command_core, project_config['freebayes_vcf'], bam_file_list))

        project_config['working_vcf'] = project_config['freebayes_vcf']
        instructions.append((command, logfile))
    else:
        sys.stderr.write("ERROR: Mode: %s not supported\n" % project_config['mode'])
        sys.exit()

    sys.stdout.write("Running FreeBayes for samples in project %s\n" % project_config['project_name'])
    pipe.execute_multiprocess(instructions, int(tool_config['freebayes']['num_cores']))
    sys.stdout.write("Finished Running FreeBayes\n")


def run_freebayes_pooled():
    """Use a pooled normal sample with freebayes"""

    raise NotImplementedError


def run_freebayes_matched():
    """Run FreeBayes in matched tumor/normal mode"""

    raise NotImplementedError()
