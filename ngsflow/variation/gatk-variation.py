"""
.. module:: GATK-Variation
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling GATK UnifiedGenotyper and HaplotypeCaller.

.. moduleauthor:: Daniel Gaston <daniel.gaston@gmail.com>


"""

__author__ = 'dgaston'

import os
import sys

from helenus.workflow import pipeline as pipe
from helenus.workflow.elements.utilities import bgzip_and_tabix_vcf_instructions, bgzip_and_tabix_vcf


def run_haplotypecaller(project_config, sample_config, tool_config, resource_config):
    """Call and Genotype variants using the HaplotypeCaller"""

    instructions = list()
    gvcfs = list()
    project_config['haplotypecaller_vcf'] = "%s.hc.vcf" % project_config['project_name']
    if not project_config['ensemble_calling']:
        project_config['cohort_vcf'] = project_config['haplotypecaller_vcf']

    for sample in sample_config:
        sample['gvcf'] = "%s.gvcf" % sample['name']
        logfile = "%s.hc.log" % sample['name']

        if os.path.isfile(sample['gvcf']):
            sys.stderr.write("WARNING: File %s already exists, defaulting to not overwriting. Halt execution and "
                             "delete file if regeneration is required" % sample['gvcf'])
        else:
            command = ("java -Xmx%sg -jar %s -T HaplotypeCaller -R %s -I %s -o %s --emitRefConfidence GVCF "
                       "--variant_index_type LINEAR --variant_index_parameter 128000 --dbsnp %s" %
                       (tool_config['gatk-haplotypecaller']['max_mem'], tool_config['gatk']['bin'],
                        resource_config['reference_genome'], sample['working_bam'], sample['gvcf'],
                        resource_config['dbsnp']))

            instructions.append((command, logfile))
        gvcfs.append("--variant %s" % sample['gvcf'])

    sys.stdout.write("Running HaplotypeCaller for samples in project %s\n" % project_config['project_name'])
    pipe.execute_multiprocess(instructions, int(tool_config['gatk-haplotypecaller']['num_cores']))
    sys.stdout.write("Finished Running HaplotypeCaller\n")

    sys.stdout.write("Running Joint Genotyping on Cohort\n")
    gvcf_string = " ".join(gvcfs)
    logfile = "%s.jointgenotyping.log" % project_config['project_name']
    command = ("java -Xmx%sg -jar %s -T GenotypeGVCFs -R %s %s -o %s" %
               (tool_config['gatk-haplotypecaller']['max_mem'], tool_config['gatk']['bin'],
                resource_config['reference_genome'], gvcf_string, resource_config['haplotypecaller_vcf']))
    code = pipe.run_and_log_command(command, logfile)
    pipe.check_return_codes(code)

    project_config['working_vcf'] = project_config['haplotypecaller_vcf']

    sys.stdout.write("Finished Joint Genotyping Cohort\n")


def run_unifiedgenotyper(project_config, sample_config, tool_config, resource_config):
    """Call Multi-Sample Variants with UnifiedGenotyper"""

    sample_inputs = list()

    for sample in sample_config:
        sample = "-I %s" % sample['working_bam']
        sample_inputs.append(sample)
    sample_bam_string = " ".join(sample_inputs)

    sys.stdout.write("Running UnifiedGenotyper for samples in project %s\n" % project_config['project_name'])
    project_config['unifiedgenotyper'] = "%s.ug.vcf" % project_config['project_name']
    if not project_config['ensemble_calling']:
        project_config['cohort_vcf'] = project_config['unifiedgenotyper_vcf']

    logfile = "%s.unifiedgenotyper.log" % project_config['project_name']

    command = ("java -Xmx%sg -jar %s -T UnifiedGenotyper -R %s %s -o %s --dbsnp %s -stand_call_conf 50.0 "
               "-stand_emit_conf 10.0 -dcov 200 --genotype_likelihoods_model BOTH --output_mode EMIT_VARIANTS_ONLY "
               "-nt %s" % (tool_config['gatk-unifiedgenotyper']['max_mem'], tool_config['gatk']['bin'],
                           resource_config['reference_genome'], sample_bam_string,
                           project_config['unifiedgenotyper_vcf'], resource_config['dbsnp'],
                           tool_config['gatk-unifiedgenotyper']['num_cores']))

    code = pipe.run_and_log_command(command, logfile)
    pipe.check_return_code(code)

    project_config['working_vcf'] = project_config['unifiedgenotyper_vcf']
    sys.stdout.write("Finished UnifiedGenotyper\n")

