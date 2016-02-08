__author__ = 'dgaston'

from helenus.workflow import pipeline as pipe
from helenus.workflow.elements.utilities import bgzip_and_tabix_vcf_instructions, bgzip_and_tabix_vcf


def run_fermikit(project_config, sample_config, tool_config, resource_config):
    """Run fermi kit assembly-based mapping and variant calling on samples"""

    for sample in sample_config:
        sample['fermi_make'] = "%s.mak" % sample['name']
        setup_log = "%s.fermi_setup.log" % sample['name']
        make_log = "%s.fermi_make.log" % sample['name']
        run_log = "%s.fermi_calling.log" % sample['name']

        setup_command = ("%s unitig -s%s -l%s -p %s -t%s \"seqtk mergepe "
                         "%s %s | trimadap-mt -p%s\" > %s" %
                         (tool_config['fermi']['fermi2.pl'], project_config['genome_size'],
                          project_config['read_size'], sample['name'], tool_config['fermi']['num_cores'],
                          sample['fastq1'], sample['fastq2'], tool_config['fermi']['num_cores'], sample['fermi_make']))

        make_command = ("make -f %s" % (sample['fermi_make']))

        calling_command = ("%s -t%s %s %s.mag.gz | sh" % (tool_config['fermi']['run-calling'],
                                                          tool_config['fermi']['num_cores'],
                                                          resource_config['reference_genome'], sample['name']))

        code = pipe.run_and_log_command(setup_command, setup_log)
        pipe.check_return_code(code)

        code = pipe.run_and_log_command(make_command, make_log)
        pipe.check_return_code(code)

        code = pipe.run_and_log_command(calling_command, run_log)
        pipe.check_return_code(code)
