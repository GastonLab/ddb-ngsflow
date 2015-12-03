"""
.. module:: annotation
   :platform: Unix, OSX
   :synopsis: A module of methods and utilities for annotating files, post-processing annotated VCFs, and transforming
   into additional formats.

.. moduleauthor:: Daniel Gaston <daniel.gaston@gmail.com>


"""

from ngsflow import pipeline


def snpeff(job, config, sample, input_vcf):
    """Annotate the specified VCF using snpEff
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_vcf: The input_vcf file name to process.
    :type input_vcf: str.
    :returns:  str -- The output vcf file name.
    """

    output_vcf = "{}.snpEff.{}.vcf".format(sample, config['snpeff']['reference'])
    logfile = "{}.snpeff.log".format(sample)

    snpeff_command = ("java",
                      "-Xmx{}g".format(config['snpeff']['max_mem']),
                      "-jar",
                      "{}".format(config['snpeff']['bin']),
                      "-v",
                      "{}".format(config['snpeff']['reference']),
                      "{}".format(input_vcf),
                      ">"
                      "{}".format(output_vcf))

    job.fileStore.logToMaster("snpEff Command: {}\n".format(snpeff_command))
    pipeline.run_and_log_command(" ".join(snpeff_command), logfile)

    return output_vcf


def gemini(job, config, sample, input_vcf):
    """Take the specified VCF and use GEMINI to add additional annotations and convert to database format
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_vcf: The input_vcf file name to process.
    :type input_vcf: str.
    :returns:  str -- The output GEMINI database name.
    """

    db = "{}.snpEff.{}.db".format(sample, config['snpeff']['reference'])
    logfile = "{}.gemini.log".format(sample)

    command = ("{}".format(config['gemini']['bin']),
               "load",
               "--cores",
               "{}".format(config['gemini']['num_cores']),
               "-v",
               "{}".format(input_vcf),
               "-t",
               "snpEff",
               "{}".format(db))

    job.fileStore.logToMaster("GEMINI Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return db
