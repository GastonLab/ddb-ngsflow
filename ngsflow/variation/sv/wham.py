__author__ = 'dgaston'

import sys

from helenus.workflow import pipeline as pipe
from helenus.workflow.elements.utilities import bgzip_and_tabix_vcf_instructions, bgzip_and_tabix_vcf


def run_wham(project_config, sample_config, tool_config, resource_config):
    """Run WHAM on a single sample"""

    instructions = list()
    instructions2 = list()

    for sample in sample_config:
        logfile1 = "%s.wham.raw.log" % sample['name']
        errfile1 = "%s.wham.raw.err"

        logfile2 = "%s.wham.class.log" % sample['name']
        errfile2 = "%s.wham.class.err"

        logfile3 = "%s.wham.sort.log" % sample['name']

        sample['wham_raw'] = "%s.wham.raw.vcf" % sample['name']
        sample['class_vcf'] = "%s.wham.class.vcf" % sample['name']
        sample['wham_vcf'] = "%s.sorted.wham.class.vcf" % sample['name']

        command1 = ("%s -f %s -e %s -t %s > %s 2> %s" % (tool_config['wham']['bin'],
                                                         resource_config['reference_genome'],
                                                         resource_config['merged_regions'], sample['working_bam'],
                                                         sample['wham_raw'], errfile1))

        command2 = ("python %s %s %s > %s 2> %s" % (tool_config['wham']['class_bin'], sample['wham_raw'],
                                                    tool_config['wham']['wham_training'], sample['class_vcf'],
                                                    errfile2))

        command3 = ("cat %s | %s -c > %s" % (sample['class_vcf'], tool_config['vcftools_sort']['bin'],
                                             sample['wham_vcf']))

        instructions.append((command2, logfile2))
        instructions2.append((command3, logfile3))

        sys.stdout.write("Running WHAM on sample\n")
        code = pipe.run_and_log_command(command1, logfile1)
        pipe.check_return_code(code)

    sys.stdout.write("Running WHAM\n")
    pipe.execute_multiprocess(instructions, int(tool_config['wham']['num_cores']))

    sys.stdout.write("Sorting WHAM VCF files\n")
    pipe.execute_multiprocess(instructions2, int(tool_config['wham']['num_cores']))
    sys.stdout.write("Finished WHAM\n")
