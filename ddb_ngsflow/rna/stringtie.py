"""
.. module:: stringtie
   :platform: Unix, OSX
   :synopsis: A module of methods for working with the StringTie transcript assembly and quantification program
   into additional formats.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>


"""

from ddb_ngsflow import pipeline


def add_additional_options(command_list, config, flags):
    if 'keep_retained' in flags:
        command_list.append("-i")

    return command_list


def stringtie_first(job, config, name, samples, flags):
    """Perform transcript assembly and quantification with StringTie
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param samples: The samples info and config dictionary.
    :type samples: dict.
    :param flags: Flags for extra parameter settings.
    :type flags: list.
    :returns:  str -- The transcript assembly GTF file name.
    """

    logfile = "{}.stringtie_first.log".format(name)
    outfile = "{}.stringtie_first.gtf".format(name)

    command = ["{}".format(config['stringtie']['bin']),
               "{}".format(samples[name]['bam']),
               "-p {}".format(config['stringtie']['num_cores']),
               "-G {}".format(config['transcript_reference_gff']),
               "-f 0.05",
               "-m 100",
               "-o {}".format(outfile)
               ]

    command = add_additional_options(command, config, flags)

    job.fileStore.logToMaster("StringTie Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return outfile


def stringtie(job, config, name, samples, flags):
    """Perform transcript assembly and quantification with StringTie
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param samples: The samples info and config dictionary.
    :type samples: dict.
    :param flags: Flags for extra parameter settings.
    :type flags: list.
    :returns:  str -- The transcript assembly GTF file name.
    """

    logfile = "{}.stringtie.log".format(name)
    outfile = "{}.stringtie.gtf".format(name)
    abundances_file = "{}.gene_abundances.txt".format(name)

    command = ["{}".format(config['stringtie']['bin']),
               "{}".format(samples[name]['bam']),
               "-p {}".format(config['stringtie']['num_cores']),
               "-G {}".format(config['merged_transcript_reference']),
               "-A {}".format(abundances_file),
               "-f 0.05",
               "-m 100",
               "-B",
               "-e",
               "-o {}".format(outfile)
               ]

    command = add_additional_options(command, config, flags)

    job.fileStore.logToMaster("StringTie Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return outfile


def stringtie_merge(job, config, samples, flags, transcripts_list):
    """Perform transcript assembly and quantification with StringTie
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param samples: The samples info and config dictionary.
    :type samples: dict.
    :param flags: Flags for extra parameter settings.
    :type flags: list.
    :returns:  str -- The transcript assembly GTF file name.
    """

    logfile = "{}.stringtie_merge.log".format(config['run_id'])
    outfile = "{}.stringtie.merged.gtf".format(config['run_id'])

    command = ["{}".format(config['stringtie']['bin']),
               "{}".format(transcripts_list),
               "--merge",
               "-p {}".format(config['stringtie']['num_cores']),
               "-G {}".format(config['merged_transcripts']),
               "-o {}".format(outfile)
               ]

    command = add_additional_options(command, config, flags)

    job.fileStore.logToMaster("StringTie Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return outfile
