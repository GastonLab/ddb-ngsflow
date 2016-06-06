"""
.. module:: Pisces (Illumina Tumor-Only Somatic Variant Caller)
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling Pisces.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>
"""

from ddb_ngsflow import pipeline


def pisces(job, config, name, input_bam):
    """Run Pisces on a single sample
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    output_vcf = "{}.pisces.vcf".format(name)
    logfile = "{}.pisces.log".format(name)
    command = ["{}".format(config['pisces']['bin']),
               "-B",
               "-t",
               "{}".format(config['pisces']['num_cores']),
               "-ThreadByChr",
               "{}".format(input_bam),
               "-g",
               "{}".format(config['reference']),
               "-f",
               "{}".format(config['min_alt_af']),
               "-b",
               "{}".format(config['min_bq']),
               "-fo",
               "False",
               "-q",
               "{}".format(config['max_var_qscore']),
               "-c",
               "{}".format(config['coverage_threshold']),
               "-s",
               "{}".format(config['sb_threshold']),
               "-a",
               "{}".format(config['min_var_qscore']),
               "-F",
               "{}".format(config['var_qscore_threshold']),
               "-gVCF",
               "True"]

    job.fileStore.logToMaster("Pisces Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return output_vcf
