__author__ = 'dgaston'

import multiprocessing

from ngsflow import pipeline


def snpeff(job, config, sample, input_vcf, max_mem):
    """Run snpEff Annotations"""

    output_vcf = "{}.snpEff.{}.vcf".format(sample, config['snpeff']['reference'])
    logfile = "{}.snpeff.log".format(sample)

    snpeff_command = ("java",
                      "-Xmx{}g".format(max_mem),
                      "-jar",
                      "{}".format(config['snpeff']['bin']),
                      "-classic",
                      "-formatEff",
                      "-v",
                      "{}".format(config['snpeff']['reference']),
                      "{}".format(input_vcf),
                      "{}".format(output_vcf))

    job.fileStore.logToMaster("snpEff Command: {}\n".format(snpeff_command))
    pipeline.run_and_log_command(" ".join(snpeff_command), logfile)

    return output_vcf


def gemini(job, config, sample, input_vcf):
    """Run GEMINI on a per sample basis"""

    db = "{}.snpEff.{}.db".format(sample, config['snpeff']['reference'])
    logfile = "{}.gemini.log".format(sample)

    command = ("{}".format(config['gemini']),
               "load",
               "--cores",
               "{}".format(multiprocessing.cpu_count()),
               "-v",
               "{}".format(input_vcf),
               "-t",
               "snpEff",
               "{}".format(db))

    job.fileStore.logToMaster("GEMINI Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return db


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
