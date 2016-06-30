"""
.. module:: hisat
   :platform: Unix, OSX
   :synopsis: A module of methods for working with the HiSat2 RNA alignment program
   into additional formats.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>


"""

import os
from ddb_ngsflow import pipeline


def add_additional_options(command_list, config, flags):
    if 'stranded' in flags:
        command_list.append("--rna-strandness {}".format(config['library-type']))

    if 'max_intron' in flags:
        command_list.append("--max-intronlen {}".format(config['hisat']['max_intron_size']))

    return command_list


def hisat_paired(job, config, name, samples, flags):
    """Align RNA-Seq data to a reference using HiSat2
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param samples: The samples info and config dictionary.
    :type samples: dict.
    :returns:  str -- The output bam file name.
    """

    logfile = "{}.hisat.log".format(name)
    output = "{}.hisat.sorted.bam".format(name)
    temp = "{}.hisat.sort.temp".format(name)

    hisat_cmd = ["{}".format(config['hisat']['bin']),
                 "-p {}".format(config['hisat']['num_cores']),
                 "--dta",
                 "-x {}".format(config['hisat']['index']),
                 "-1 {}".format(samples[name]['fastq1']),
                 "-2 {}".format(samples[name]['fastq2'])
                 ]

    hisat_cmd = add_additional_options(hisat_cmd, config, flags)

    view_cmd = ["{}".format(config['samtools']['bin']),
                "view",
                "-u",
                "-"]

    sort_cmd = ["{}".format(config['samtools']['bin']),
                "sort",
                "-@",
                "{}".format(config['hisat']['num_cores']),
                "-O",
                "bam",
                "-o",
                "{}".format(output),
                "-T",
                "{}".format(temp),
                "-"]

    command = "{} | {} | {}".format(" ".join(hisat_cmd), " ".join(view_cmd), " ".join(sort_cmd))

    job.fileStore.logToMaster("HiSat2 Command: {}\n".format(command))
    pipeline.run_and_log_command(command, logfile)

    return output


def hisat_unpaired(job, config, name, samples, flags):
    """Align RNA-Seq data to a reference using HiSat2
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param samples: The samples info and config dictionary.
    :type samples: dict.
    :returns:  str -- The output bam file name.
    """

    working_dir = os.getcwd()

    logfile = "{}.hisat.log".format(name)
    output = "{}.hisat.sorted.bam".format(name)
    unaligned = os.path.join(working_dir, "{}.unaligned.sam".format(name))
    temp = "{}.hisat.sort.temp".format(name)

    hisat_cmd = ["{}".format(config['hisat']['bin']),
                 "-p {}".format(config['hisat']['num_cores']),
                 "--dta",
                 "-x {}".format(config['hisat']['index']),
                 "-U {}".format(samples[name]['fastq1']),
                 "--un {}".format(unaligned)
                 ]

    view_cmd = ["{}".format(config['samtools']['bin']),
                "view",
                "-u",
                "-"]

    sort_cmd = ["{}".format(config['samtools']['bin']),
                "sort",
                "-@",
                "{}".format(config['hisat']['num_cores']),
                "-O",
                "bam",
                "-o",
                "{}".format(output),
                "-T",
                "{}".format(temp),
                "-"]

    command = "{} | {} | {}".format(" ".join(hisat_cmd), " ".join(view_cmd), " ".join(sort_cmd))

    job.fileStore.logToMaster("HiSat2 Command: {}\n".format(command))
    pipeline.run_and_log_command(command, logfile)

    return output
