"""
.. module:: LoFreq
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling LoFreq.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>

"""

from ddb_ngsflow import pipeline


def run_lowfreq(job, config, sample, input_bam):
    """Run LoFreq on an an unmatched tumour sample and call somatic variants
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    vcf = "{}.lofreq.vcf".format(sample)
    logfile = "{}.lofreq.log".format(sample)

    command = ("{}".format(config['lofreq']['bin']),
               "somatic",
               "-t",
               "{}".format(input_bam),
               "--call-indels"
               "-f",
               "{}".format(config['reference']),
               "--threads",
               "{}".format(config['lofreq']['num_cores']),
               "-d",
               "{}".format(config['dbsnp']),
               "-o",
               "{}".format(vcf))

    job.fileStore.logToMaster("LoFreq Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)
