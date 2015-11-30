"""
.. module:: indelminer
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling indelminer.

.. moduleauthor:: Daniel Gaston <daniel.gaston@gmail.com>


"""

__author__ = 'dgaston'

from ngsflow import pipeline


# indeliminer can be quite slow, doesn't appear to have inbuilt multi-threading, and can use a lot of RAM (>4GB)
def indelminer_single(job, config, sample, input_bam):
    """Run indelminer for calling indels in sample data"""

    indelminer_vcf = "{}.indelminer.vcf".format(sample)
    logfile = "{}.indelminer.log".format(sample)

    command = ("{}".format(config['indelminer']['bin']),
               "{}".format(config['reference']),
               "{}".format(input_bam),
               ">",
               "{}".format(indelminer_vcf))

    job.fileStore.logToMaster("IndelMiner Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return indelminer_vcf
