"""
.. module:: platypus
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling platypus.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>

"""

from ddb_ngsflow import pipeline


# This needs to be fixed, new regions need to be defined for all targeted panels to use this
def platypus_single(job, config, name, samples, input_bam):
    """Run Platypus on an an unmatched tumour sample and call somatic variants
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    platypus_vcf = "{}.platypus.vcf".format(name)
    platypus_log = "{}.platypus.log".format(name)
    internal_log = "{}.platypus_internal.log".format(name)

    platypus_command = ("{}".format(config['platypus']['bin']),
                        "callVariants",
                        "--refFile={}".format(config['reference']),
                        "--regions={}".format(samples[name]['regions']),
                        "--assemble=1",
                        "--assembleBadReads=1",
                        "--assembleBrokenPairs=1",
                        "--filterDuplicates=0",
                        "--minVarFreq={}".format(config['min_alt_af']),
                        "--nCPU={}".format(config['platypus']['num_cores']),
                        "--logFileName={}".format(internal_log),
                        "--bamFiles={}".format(input_bam),
                        "--output={}".format(platypus_vcf))

    job.fileStore.logToMaster("Platypus Command: {}\n".format(platypus_command))
    pipeline.run_and_log_command(" ".join(platypus_command), platypus_log)

    return platypus_vcf
