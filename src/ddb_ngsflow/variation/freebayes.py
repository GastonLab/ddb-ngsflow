"""
.. module:: freebayes
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling FreeBayes.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>
"""

from ngsflow import pipeline


def freebayes_single(job, config, sample, input_bam):
    """Run FreeBayes without a matched normal sample
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    freebayes_vcf = "{}.freebayes.vcf".format(sample)
    logfile = "{}.freebayes.log".format(sample)

    command = ("{}".format(config['freebayes']['bin']),
               "--fasta-reference",
               "{}".format(config['reference']),
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
    pipeline.run_and_log_command(" ".join(command), logfile)

    return freebayes_vcf


def freebayes_pooled():
    """Use a pooled normal sample with freebayes"""

    raise NotImplementedError


def freebayes_matched():
    """Run FreeBayes in matched tumor/normal mode"""

    raise NotImplementedError()
