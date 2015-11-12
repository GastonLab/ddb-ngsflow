__author__ = 'dgaston'

import time

from ngsflow import pipeline


def run_freebayes_single(job, config, sample, input_bam):
    """Run FreeBayes without a matched normal sample"""

    freebayes_vcf = "{}.freebayes.vcf".format(sample)
    logfile = "{}.freebayes.log".format(sample)

    command = ("{}".format(config['freebayes']),
               "--fasta-reference",
               "{}".format(config['reference']),
               "-t",
               "{}".format(config['merged_regions']),
               "--min-alternate-fraction",
               "{}".format(config['min_alt_af']),
               "--pooled-discrete",
               "--pooled-continuous",
               "--genotype-qualities",
               "--report-genotype-likelihood-max",
               "--allele-balance-priors-off",
               "--use-duplicate-reads",
               "--min-repeat-entropy 1",
               "-v",
               "{}".format(freebayes_vcf),
               "{}".format(input_bam))

    job.fileStore.logToMaster("FreeBayes Command: {}\n".format(command))
    # pipeline.run_and_log_command(" ".join(command), logfile)
    time.sleep(2)
    return freebayes_vcf


def run_freebayes_pooled():
    """Use a pooled normal sample with freebayes"""

    raise NotImplementedError


def run_freebayes_matched():
    """Run FreeBayes in matched tumor/normal mode"""

    raise NotImplementedError()
