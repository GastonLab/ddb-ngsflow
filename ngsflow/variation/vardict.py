__author__ = 'dgaston'

import time
import multiprocessing

from ngsflow import pipeline
from ngsflow.utils import utilities


def vardict_matched():
    """Run VarDict in matched tumor/normal mode"""

    raise NotImplementedError()


def vardict_single(job, config, sample, input_bam):
    """Run VarDict without a matched normal sample"""

    vardict_vcf = "{}.vardict.vcf".format(sample)
    logfile = "{}.vardict.log".format(sample)

    vardict = ("{}".format(config['vardict']['bin']),
               "-G",
               "{}".format(config['reference']),
               "-z",
               "-c",
               "1",
               "-S",
               "2",
               "-E",
               "3",
               "-g",
               "4",
               "-th",
               "{}".format(multiprocessing.cpu_count()),
               "-f",
               "{}".format(config['min_alt_af']),
               "-N",
               "{}".format(sample),
               "-b",
               "{}".format(input_bam),
               "{}".format(config['regions']))

    vardict2vcf = ("{}".format(config['vardict2vcf']),
                   "-E",
                   "-f",
                   "{}".format(config['min_alt_af']),
                   "-N",
                   "{}".format(sample))

    command = ("{vardict} | {strandbias} | {vardict2vcf} > {vcf}".format(vardict=vardict,
                                                                         strandbias=config['vardict_strandbias']['bin'],
                                                                         vardict2vcf=config['vardict2vcf']['bin'],
                                                                         vcf=vardict_vcf))

    job.fileStore.logToMaster("VarDict Command: {}\n".format(command))
    # pipeline.run_and_log_command(command, logfile)
    utilities.touch(vardict_vcf)
    time.sleep(2)

    return vardict_vcf
