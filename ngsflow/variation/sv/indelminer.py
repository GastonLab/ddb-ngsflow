__author__ = 'dgaston'

import time

from ngsflow import pipeline


def run_indelminer_single(job, config, sample, input_bam):
    """Run indelminer for calling indels in sample data"""

    indelminer_vcf = "{}.indelminer.vcf".format(sample)
    logfile = "{}.indelminer.log".format(sample)

    command = ("{}".format(config['indelminer']),
               "{}".format(config['reference']),
               "{}".format(input_bam),
               ">",
               "{}".format(indelminer_vcf))

    job.fileStore.logToMaster("IndelMiner Command: {}\n".format(command))
    # pipeline.run_and_log_command(" ".join(command), logfile)
    time.sleep(2)

    return indelminer_vcf
