"""
.. module:: star
   :platform: Unix, OSX
   :synopsis: A module of methods for working with the STAR RNA alignment program
   into additional formats.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>


"""

from ddb_ngsflow import pipeline


def add_additional_options(command_list, config, flags):
    if 'compressed' in flags:
        command_list.append("--readFilesCommand {}".format(config['compression']))

    if 'encode_options' in flags:
        encode_options = ("--outFilterType BySJout",
                          "--outFilterMultimapNmax 20",
                          "--alignSJDBoverhangMin 1",
                          "--outFilterMismatchNMax 999",
                          "--alignIntronMin 20",
                          "--alignIntronMax 1000000",
                          "--alignMatesGapMax 1000000")
        command_list.extend(encode_options)

    return command_list


def star_paired(job, config, name, samples, flags):
    """Align RNA-Seq data to a reference using STAR
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param samples: The samples info and config dictionary.
    :type samples: dict.
    :returns:  str -- The output vcf file name.
    """

    output = "{}.star.output".format(name)
    logfile = "{}.star.log".format(name)

    command = ("{}".format(config['star']['bin']),
               "--genomeDir {}".format(config['star']['index']),
               "--runThreadN {}".format(config['star']['num_cores']),
               "--readFilesIn {} {}".format(samples[name]['fastq1'], samples[name]['fastq2']),
               "--outFileNamePrefix {}".format(output),
               "--outReadsUnmapped Fastx"
               )

    command = add_additional_options(command, config, flags)

    job.fileStore.logToMaster("STAR Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return output


def star_unpaired(job, config, name, samples, flags):
    """Align RNA-Seq data to a reference using STAR
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param samples: The samples info and config dictionary.
    :type samples: dict.
    :returns:  str -- The output vcf file name.
    """

    output = "{}.star.output".format(name)
    logfile = "{}.star.log".format(name)

    command = ("{}".format(config['star']['bin']),
               "--genomeDir {}".format(config['star']['index']),
               "--runThreadN {}".format(config['star']['num_cores']),
               "--readFilesIn {}".format(samples[name]['fastq1']),
               "--outFileNamePrefix {}".format(output),
               "--outReadsUnmapped Fastx"
               )

    command = add_additional_options(command, config, flags)

    job.fileStore.logToMaster("STAR Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return output


# def star_paired_compressed_basic(job, config, name, samples):
#     """Align RNA-Seq data to a reference using STAR
#     :param config: The configuration dictionary.
#     :type config: dict.
#     :param name: sample name.
#     :type name: str.
#     :param samples: The samples info and config dictionary.
#     :type samples: dict.
#     :returns:  str -- The output vcf file name.
#     """
#
#     output_dir = "{}.star.output".format(name)
#     logfile = "{}.star.log".format(name)
#
#     command = ("{}".format(config['star']['bin']),
#                "--genomeDir {}".format(config['star']['index']),
#                "--runThreadN {}".format(config['star']['num_cores']),
#                "--readFilesIn {} {}".format(samples[name]['fastq1'], samples[name]['fastq2']),
#                "--readFilesCommand {}".format(config['compression']),
#                "--outFileNamePrefix {}".format(output_dir)
#                )
#
#     job.fileStore.logToMaster("STAR Command: {}\n".format(command))
#     pipeline.run_and_log_command(" ".join(command), logfile)
#
#     return output_dir
#
#
# def star_paired_uncompressed_encode(job, config, name, samples):
#     """Align RNA-Seq data to a reference using STAR
#     :param config: The configuration dictionary.
#     :type config: dict.
#     :param name: sample name.
#     :type name: str.
#     :param samples: The samples info and config dictionary.
#     :type samples: dict.
#     :returns:  str -- The output vcf file name.
#     """
#
#     output_dir = "{}.star.output".format(name)
#     logfile = "{}.star.log".format(name)
#
#     command = ("{}".format(config['star']['bin']),
#                "--genomeDir {}".format(config['star']['index']),
#                "--runThreadN {}".format(config['star']['num_cores']),
#                "--readFilesIn {} {}".format(samples[name]['fastq1'], samples[name]['fastq2']),
#                "--outFilterType BySJout",
#                "--outFilterMultimapNmax 20",
#                "--alignSJDBoverhangMin 1",
#                "--outFilterMismatchNMax 999",
#                "--alignIntronMin 20",
#                "--alignIntronMax 1000000",
#                "--alignMatesGapMax 1000000",
#                "--outFileNamePrefix {}".format(output_dir)
#                )
#
#     job.fileStore.logToMaster("STAR Command: {}\n".format(command))
#     pipeline.run_and_log_command(" ".join(command), logfile)
#
#     return output_dir
#
#
# def star_paired_compressed_encode(job, config, name, samples):
#     """Align RNA-Seq data to a reference using STAR
#     :param config: The configuration dictionary.
#     :type config: dict.
#     :param name: sample name.
#     :type name: str.
#     :param samples: The samples info and config dictionary.
#     :type samples: dict.
#     :returns:  str -- The output vcf file name.
#     """
#
#     output_dir = "{}.star.output".format(name)
#     logfile = "{}.star.log".format(name)
#
#     command = ("{}".format(config['star']['bin']),
#                "--genomeDir {}".format(config['star']['index']),
#                "--runThreadN {}".format(config['star']['num_cores']),
#                "--readFilesIn {} {}".format(samples[name]['fastq1'], samples[name]['fastq2']),
#                "--outFilterType BySJout",
#                "--outFilterMultimapNmax 20",
#                "--alignSJDBoverhangMin 1",
#                "--outFilterMismatchNMax 999",
#                "--alignIntronMin 20",
#                "--alignIntronMax 1000000",
#                "--alignMatesGapMax 1000000",
#                "--outFileNamePrefix {}".format(output_dir),
#                "--readFilesCommand {}".format(config['compression'])
#                )
#
#     job.fileStore.logToMaster("STAR Command: {}\n".format(command))
#     pipeline.run_and_log_command(" ".join(command), logfile)
#
#     return output_dir
