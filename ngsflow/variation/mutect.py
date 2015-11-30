"""
.. module:: mutect
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling MuTect.

.. moduleauthor:: Daniel Gaston <daniel.gaston@gmail.com>


"""

__author__ = 'dgaston'

import multiprocessing

from ngsflow import pipeline
from ngsflow.utils import utilities


def mutect_pon():
    """Run MuTect with a synthetic Panel of Normals"""

    raise NotImplementedError()


def mutect_matched():
    """Run MuTect on paired tumor normal data"""

    raise NotImplementedError()


def mutect_single(job, config, sample, input_bam, max_mem):
    """Run MuTect without paired normal samples. Use vcf-subset to remove none column"""

    mutect_vcf = "{}.mutect.vcf".format(sample)
    temp_mutect = "{}.tempmutect.vcf".format(sample)

    output_stats = "{}.mutectstats.txt".format(sample)
    sample_coverage = "{}.mutectcoverage.wig.txt".format(sample)

    mutect_logfile = "{}.mutect.log".format(sample)
    subset_log = "{}.mutect_subset.log".format(sample)

    mutect_command = ("java",
                      "-Xmx{}g".format(max_mem),
                      "-jar",
                      "{}".format(config['mutect']['bin']),
                      "-T",
                      "MuTect",
                      "-R",
                      "{}".format(config['reference']),
                      "--dbsnp",
                      "{}".format(config['dbsnp']),
                      "--cosmic",
                      "{}".format(config['cosmic']),
                      "--enable_extended_output",
                      "-I:tumor",
                      "{}".format(input_bam),
                      "--coverage_file",
                      "{}".format(sample_coverage),
                      "-nt",
                      "{}".format(multiprocessing.cpu_count()),
                      "-o",
                      "{}".format(output_stats),
                      "-vcf",
                      "{}".format(temp_mutect))

    subset_command = ("cat",
                      "%s".format(temp_mutect),
                      "|",
                      "{}".format(config['vcftools_subset']['bin']),
                      "-e",
                      "-c",
                      "{}".format(sample),
                      ">",
                      "{}".format(mutect_vcf))

    job.fileStore.logToMaster("MuTect Command: {}\n".format(mutect_command))
    pipeline.run_and_log_command(" ".join(mutect_command), mutect_logfile)

    job.fileStore.logToMaster("Subset Command: {}\n".format(subset_command))
    pipeline.run_and_log_command(" ".join(subset_command), subset_logfile)

    return mutect_vcf
