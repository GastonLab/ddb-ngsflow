__author__ = 'dgaston'

import os
import time
import multiprocessing

from ngsflow import pipeline
from ngsflow.utils import utilities


def scalpel_single(job, config, sample, input_bam):
    """Run scalpel for variant calling"""
    cwd = os.getcwd()
    output_dir = os.path.join(cwd, "{}-scalpel-output".format(sample))
    scalpel_vcf = os.path.join(output_dir, "variants.indel.vcf")
    logfile = "{}.scalpel.log".format(sample)

    scalpel_command = ("{}".format(config['scalpel-discovery']['bin']),
                       "--single",
                       "--intarget",
                       "--covthr",
                       "3",
                       "--lowcov",
                       "1",
                       "--ref",
                       "{}".format(config['reference']),
                       "--bed",
                       "{}".format(config['regions']),
                       "--format",
                       "vcf",
                       "--numprocs",
                       "{}".format(multiprocessing.cpu_count()),
                       "--bam",
                       "{}".format(input_bam),
                       "--dir",
                       "{}".format(output_dir))

    job.fileStore.logToMaster("Scalpel Command: {}\n".format(scalpel_command))
    # pipeline.run_and_log_command(" ".join(scalpel_command), logfile)
    utilities.touch(scalpel_vcf)
    time.sleep(2)

    return scalpel_vcf
