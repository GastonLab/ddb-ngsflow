__author__ = 'dgaston'

import sys

from .. import pipeline as pipe


def run_snpeff(project_config, sample_config, tool_config, resource_config):
    """Run snpEff Annotations"""

    instructions = list()

    command_core = ("java -Xmx%sg -jar %s -classic -formatEff -v %s" % (tool_config['snpeff']['max_mem'],
                                                                        tool_config['snpeff']['bin'],
                                                                        tool_config['snpeff']['reference']))
    if project_config['mode'] == "per_sample":
        for sample in sample_config:
            sample['snpeff_vcf'] = "%s.snpEff.%s.vcf" % (sample['name'], tool_config['snpeff']['reference'])
            logfile = "%s.snpeff.log" % sample['name']

            command = ("%s %s > %s" % (command_core, sample['working_vcf'], sample['snpeff_vcf']))
            instructions.append((command, logfile))
            sample['working_vcf'] = sample['snpeff_vcf']

    elif project_config['mode'] == "per_cohort":
        project_config['snpeff_vcf'] = "%s.snpEff.%s.vcf" % (project_config['project_name'], 
                                                             tool_config['snpeff']['reference'])
        logfile = "%s.snpeff.log" % project_config['project_name']

        command = ("%s %s > %s" % (command_core, project_config['working_vcf'], project_config['snpeff_vcf']))
        instructions.append((command, logfile))
        project_config['working_vcf'] = project_config['snpeff_vcf']
    else:
        sys.stderr.write("ERROR: Mode: %s not supported\n" % project_config['mode'])
        sys.exit()

    sys.stdout.write("Running snpEff\n")
    pipe.execute_multiprocess(instructions, int(tool_config['snpeff']['num_cores']))
    sys.stdout.write("Finished snpEff\n")


def run_gemini(project_config, sample_config, tool_config, resource_config):
    """Run GEMINI on a per sample basis"""

    if project_config['mode'] == "per_sample":
        for sample in sample_config:
            sample['gemini_db'] = "%s.snpEff.%s.db" % (sample['name'], tool_config['snpeff']['reference'])
            logfile = "%s.gemini.log" % sample['name']

            command = ("%s load --cores %s -v %s -t snpEff %s" %
                       (tool_config['gemini']['bin'], tool_config['gemini']['num_cores'], sample['working_vcf'],
                        sample['gemini_db']))

            sys.stdout.write("Running GEMINI for sample %s\n" % sample['name'])
            code = pipe.run_and_log_command(command, logfile)
            pipe.check_return_code(code)
    elif project_config['mode'] == "per_cohort":
        project_config['gemini_db'] = "%s.snpEff.%s.db" % (project_config['project_name'],
                                                           tool_config['snpeff_parameters']['reference'])
        logfile = "%s.gemini.log" % project_config['project_name']

        command = ("%s load --cores %s -v %s -t snpEff %s" %
                   (tool_config['gemini']['bin'], tool_config['gemini']['num_cores'], project_config['working_vcf'],
                    project_config['gemini_db']))

        sys.stdout.write("Running GEMINI\n")
        code = pipe.run_and_log_command(command, logfile)
        pipe.check_return_code(code)
    else:
        sys.stderr.write("ERROR: Mode: %s not supported\n" % project_config['mode'])
        sys.exit()

    sys.stdout.write("Finished GEMINI for all samples\n")


def annotate_structural_variants(project_config, sample_config, tool_config, resource_config):
    """Annotate structural variants using a provided BED file for genes, targets, exons, etc"""

    instructions = list()

    command_core = ("%s annotate_regions -b %s -t SVAR -d \"Structural Variant Annotation Regions\" "
                    % (tool_config['vt']['bin'], resource_config['regions']))

    if project_config['mode'] == "per_sample":
        for sample in sample_config:
            logfile = "%s.annotate-structvariants.log" % sample['name']
            sample['annotated_sv_vcf'] = "%s.annotated.sv.vcf" % sample['name']

            command = ("%s -o %s %s" % (command_core, sample['annotated_sv_vcf'], sample['working_sv_vcf']))

            instructions.append((command, logfile))
    elif project_config['mode'] == "per_cohort":
        logfile = "%s.annotate-structvariants.log" % project_config['project_name']
        project_config['annotated_sv_vcf'] = "%s.annotated.sv.vcf" % project_config['project_name']

        command = ("%s -o %s %s" % (command_core, project_config['annotated_sv_vcf'], project_config['working_sv_vcf']))

        instructions.append((command, logfile))
    else:
        sys.stderr.write("ERROR: Mode: %s not supported\n" % project_config['mode'])
        sys.exit()

    sys.stdout.write("Annotating structural variants\n")
    pipe.execute_multiprocess(instructions, int(tool_config['snpeff']['num_cores']))
    sys.stdout.write("Finished annotating structural variants\n")
