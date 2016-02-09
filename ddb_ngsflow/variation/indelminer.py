"""
.. module:: indelminer
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling indelminer.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>

"""

from ddb_ngsflow import pipeline


# indeliminer can be quite slow, doesn't appear to have inbuilt multi-threading, and can use a lot of RAM (>6GB?)
# No longer using as it just takes too long and is too resource hungry.
def indelminer_single(job, config, sample, input_bam):
    """Run IndelMiner on an an unmatched tumour sample and call somatic variants
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    indelminer_vcf = "{}.indelminer.vcf".format(sample)
    logfile = "{}.indelminer.log".format(sample)

    command = ("{}".format(config['indelminer']['bin']),
               "{}".format(config['reference']),
               "sample={}".format(input_bam),
               ">",
               "{}".format(indelminer_vcf))

    job.fileStore.logToMaster("IndelMiner Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return indelminer_vcf
