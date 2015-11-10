__author__ = 'dgaston'

import sys

from helenus.workflow import pipeline as pipe
from helenus.workflow.elements.utilities import bgzip_and_tabix_vcf_instructions, bgzip_and_tabix_vcf


def run_mutect_pon():
    """Run MuTect with a synthetic Panel of Normals"""

    raise NotImplementedError()


def run_mutect_matched():
    """Run MuTect on paired tumor normal data"""

    raise NotImplementedError()


def run_mutect_single(project_config, sample_config, tool_config, resource_config):
    """Run MuTect without paired normal samples. Use vcf-subset to remove none column"""

    instructions1 = list()
    instructions2 = list()

    mutect_core = ("java -Xmx%sg -jar %s -T MuTect -R %s --dbsnp %s --cosmic %s --enable_extended_output" %
                   (tool_config['mutect']['max_mem'], tool_config['mutect']['bin'],
                    resource_config['reference_genome'], resource_config['dbsnp'], resource_config['cosmic']))

    if project_config['mode'] == "per_sample":
        for sample in sample_config:
            sample['mutect_vcf'] = "%s.mutect.vcf" % sample['name']
            temp_mutect = "%s.tempmutect.vcf" % sample['name']
            if not project_config['ensemble_calling']:
                sample['raw_vcf'] = sample['mutect_vcf']
            output_stats = "%s.mutectstats.txt" % sample['name']
            sample_coverage = "%s.mutectcoverage.wig.txt" % sample['name']
            logfile = "%s.mutect.log" % sample['name']
            logfile2 = "%s.mutect_subset.log" % sample['name']

            mutect_command = ("%s -I:tumor %s --coverage_file %s -o %s -vcf %s" %
                              (mutect_core, sample['recalibrated_bam'], sample_coverage, output_stats, temp_mutect))

            subset_command = ("cat %s | %s -e -c %s > %s" % (temp_mutect, tool_config['vcftools_subset']['bin'],
                                                             sample['name'], sample['mutect_vcf']))

            instructions1.append((mutect_command, logfile))
            instructions2.append((subset_command, logfile2))
    elif project_config['mode'] == "per_cohort":
        sys.stderr.write("ERROR: Mode: %s not supported\n" % project_config['mode'])
        sys.exit()
    else:
        sys.stderr.write("ERROR: Mode: %s not supported\n" % project_config['mode'])
        sys.exit()

    sys.stdout.write("Running MuTect with tumor only\n")
    pipe.execute_multiprocess(instructions1, int(tool_config['mutect']['num_cores']))
    pipe.execute_multiprocess(instructions2, int(tool_config['mutect']['num_cores']))
    sys.stdout.write("Finished running MuTect\n")
