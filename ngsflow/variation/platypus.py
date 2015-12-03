"""
.. module:: platypus
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling platypus.

.. moduleauthor:: Daniel Gaston <daniel.gaston@gmail.com>

"""

from ngsflow import pipeline


def platypus_single(job, config, sample, input_bam):
    """Run Platypus on an an unmatched tumour sample and call somatic variants
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    platypus_vcf = "{}.platypus.vcf".format(sample)
    platypus_log = "{}.platypus.log".format(sample)
    internal_log = "{}.platypus_internal.log".format(sample)

    platypus_command = ("{}".format(config['platypus']['bin']),
                        "callVariants",
                        "--refFile={}".format(config['reference']),
                        "--regions={}".format(config['platypus']['regions']),
                        "--assemble=1",
                        "--assembleBrokenPairs=1",
                        "--filterDuplicates=0",
                        "--nCPU={}".format(config['platypus']['num_cores']),
                        "--logFileName={}".format(internal_log),
                        "--bamFiles={}".format(input_bam),
                        "--output={}".format(platypus_vcf))

    job.fileStore.logToMaster("Platypus Command: {}\n".format(platypus_command))
    pipeline.run_and_log_command(" ".join(platypus_command), platypus_log)

    return platypus_vcf
