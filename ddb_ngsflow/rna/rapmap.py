"""
.. module:: star
   :platform: Unix, OSX
   :synopsis: A module of methods for working with the STAR RNA alignment program
   into additional formats.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>


"""

from ddb_ngsflow import pipeline
from exceptions import NotImplementedError


def rapmap_quasi_unpaired(job, config, name, samples, flags):
    """Run RapMap Quasi-Mapping procedure on unpaired sequencing data
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param samples: The samples info and config dictionary.
    :type samples: dict.
    :returns:  str -- The output vcf file name.
    """

    output = "{}".format(name)
    logfile = "{}".format(name)

    command = ["{} quasimap".format(config['rapmap']['bin']),
               "-t {}".format(config['rapmap']['num_cores']),
               "-i {}".format(config['rapmap']['index']),
               "-r {}".format(samples[name]['fastq1']),
               "-o {}".format(output)
               ]

    job.fileStore.logToMaster("RapMap Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return output


def rapmap_quasi_paired(job, config, name, samples, flags):
    """Run RapMap Quasi-Mapping procedure on paired-end sequencing data
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param samples: The samples info and config dictionary.
    :type samples: dict.
    :returns:  str -- The output vcf file name.
    """

    output = "{}".format(name)
    logfile = "{}".format(name)

    command = ["{} quasimap".format(config['rapmap']['bin']),
               "-t {}".format(config['rapmap']['num_cores']),
               "-i {}".format(config['rapmap']['index']),
               "-1 {}".format(samples[name]['fastq1']),
               "-2 {}".format(samples[name]['fastq2']),
               "-o {}".format(output)
               ]

    job.fileStore.logToMaster("RapMap Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return output

    return NotImplementedError


def rapmap_pseudo_unpaired(job, config, name, samples, flags):
    """Run RapMap Pseudo-Mapping procedure on unpaired sequencing data
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param samples: The samples info and config dictionary.
    :type samples: dict.
    :returns:  str -- The output vcf file name.
    """

    output = "{}".format(name)
    logfile = "{}".format(name)

    command = []

    # job.fileStore.logToMaster("RapMap Command: {}\n".format(command))
    # pipeline.run_and_log_command(" ".join(command), logfile)

    return NotImplementedError


def rapmap_pseudo_paired(job, config, name, samples, flags):
    """Run RapMap Pseudo-Mapping procedure on paired-end sequencing data
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param samples: The samples info and config dictionary.
    :type samples: dict.
    :returns:  str -- The output vcf file name.
    """

    output = "{}".format(name)
    logfile = "{}".format(name)

    command = []

    # job.fileStore.logToMaster("RapMap Command: {}\n".format(command))
    # pipeline.run_and_log_command(" ".join(command), logfile)

    return NotImplementedError
