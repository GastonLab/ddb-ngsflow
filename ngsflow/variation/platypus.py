__author__ = 'dgaston'

import multiprocessing

from ngsflow import pipeline


def platypus_single(job, config, sample, input_bam):
    """Run platypus on a single sample"""

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
                        "--nCPU={}".format(multiprocessing.cpu_count()),
                        "--logFileName={}".format(internal_log),
                        "--bamFiles={}".format(input_bam),
                        "--output={}".format(platypus_vcf))

    job.fileStore.logToMaster("Platypus Command: {}\n".format(platypus_command))
    pipeline.run_and_log_command(" ".join(platypus_command), platypus_log)

    return platypus_vcf
