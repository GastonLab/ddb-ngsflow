"""
.. module:: freebayes
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling ScanIndel.

.. moduleauthor:: Daniel Gaston <daniel.gaston@gmail.com>
"""

from ngsflow import pipeline


def scanindel(job, config, sample, samples, input_bam):
    """Run MANTA caller for Structural Variant Detection
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    output_vcf = "{}.merged.indel.vcf".format(sample)
    logfile = "{}.scanindel.log".format(sample)
    sample_config_file = "{}.scanindel_sample_config.txt".format(sample)

    with open(sample_config_file, 'w') as sample_config:
        sample_config.write("{id}\t{file}".format(id=sample, file=input_bam))

    command = ("{}".format(config['scanindel']['bin']),
               "-i",
               "{}".format(sample_config_file),
               "-p",
               "{}".format(config['scanindel']['config_file']),
               "--bam",
               "-F",
               "{}".format(config['min_alt_af']),
               "-t",
               "{}".format(samples[sample]['regions']))

    job.fileStore.logToMaster("Manta Configuration Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return output_vcf
