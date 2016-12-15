"""
.. module:: utilities
   :platform: Unix, OSX
   :synopsis: A module of methods for various sambamba functions. This is a general catch all module for methoods
   not better implemented elsewhere.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>

"""

import sys
import csv
from collections import defaultdict
from ddb_ngsflow import pipeline


def sambamba_region_coverage(job, config, name, samples, input_bam):
    """Run SamBambam to calculate the coverage of targeted regions
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample/library name.
    :type name: str.
    :param input_bam: The input_bam file name to process.
    :type samples: dict
    :param samples: The samples configuration dictionary
    :type input_bam: str.
    :returns:  str -- The output BED file name.
    """

    output = "{}.sambamba_coverage.bed".format(name)
    logfile = "{}.sambamba_coverage.log".format(name)

    command = ["{}".format(config['sambamba']['bin']),
               "depth region",
               "-L",
               "{}".format(samples[name]['regions']),
               "-t",
               "{}".format(config['sambamba']['num_cores']),
               "-T",
               "{}".format(config['coverage_threshold']),
               "-T",
               "{}".format(config['coverage_threshold2']),
               "{}".format(input_bam),
               ">",
               "{}".format(output)]

    job.fileStore.logToMaster("SamBamba Coverage Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return output
