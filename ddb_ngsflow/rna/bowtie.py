"""
.. module:: cufflinks
   :platform: Unix, OSX
   :synopsis: A module of methods for working with the bowtie alignment program
   into additional formats.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>


"""

from ddb_ngsflow import pipeline


def add_additional_options(command_list, config, flags):
    if 'local' in flags:
        command_list.append("--local")

    return command_list


def bowtie_unpaired(job, config, name, samples, flags):
    """Align RNA-Seq data to a reference using Bowtie2
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param samples: The samples info and config dictionary.
    :type samples: dict.
    :returns:  str -- The output vcf file name.
    """

    output = "{}.bowtie.sam".format(name)
    logfile = "{}.bowtie.log".format(name)

    if "2-stage" in flags:
        samples[name]['fastq1'] = samples[name]['unmapped_fastq']

    command = ["{}".format(config['bowtie']['bin']),
               "-x {}".format(config['bowtie']['index']),
               "-p {}".format(config['bowtie']['num_cores']),
               "-U {}".format(samples[name]['fastq1']),
               "-S {}".format(output)
               ]

    command = add_additional_options(command, config, flags)

    job.fileStore.logToMaster("Bowtie Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return output


def bowtie_paired(job, config, name, samples, flags):
    """Align RNA-Seq data to a reference using Bowtie2
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param samples: The samples info and config dictionary.
    :type samples: dict.
    :returns:  str -- The output vcf file name.
    """

    output = "{}.bowtie.sam".format(name)
    logfile = "{}.bowtie.log".format(name)

    command = ["{}".format(config['bowtie']['bin']),
               "-x {}".format(config['bowtie']['index']),
               "-p {}".format(config['bowtie']['num_cores']),
               "-1 {}".format(samples[name]['fastq1']),
               "-2 {}".format(samples[name]['fastq2']),
               "-S {}".format(output)
               ]

    command = add_additional_options(command, config, flags)

    job.fileStore.logToMaster("Bowtie Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return output
