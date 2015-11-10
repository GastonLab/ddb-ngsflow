__author__ = 'dgaston'

import os
import sys

from helenus.workflow import pipeline as pipe


def run_scalpel_single(project_config, sample_config, tool_config, resource_config):
    """Run scalpel for variant calling"""

    scalpel_command_core = ("%s --single --intarget --covthr 3 --lowcov 1 --ref %s --bed %s --format vcf "
                            "--numprocs %s" %
                            (tool_config['scalpel-discovery']['bin'], resource_config['reference_genome'],
                             resource_config['merged_regions'], tool_config['scalpel']['num_cores']))

    for sample in sample_config:
        cwd = os.getcwd()
        output_dir = "%s/%s-scalpel-output" % (cwd, sample['name'])
        logfile = "%s.scalpel.log" % sample['name']
        sample['scalpel_vcf'] = "%s/variants.indel.vcf" % output_dir

        if not project_config['ensemble_calling']:
            sample['raw_vcf'] = sample['scalpel_vcf']

        sys.stdout.write("Running Scalpel for sample %s\n" % sample['name'])
        command = ("%s --bam %s --dir %s" % (scalpel_command_core, sample['working_bam'], output_dir))
        code = pipe.run_and_log_command(command, logfile)
        pipe.check_return_code(code)
        sys.stdout.write("Finished Scalpel for sample %s\n" % sample['name'])
