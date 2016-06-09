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
        encode_options = ["--outFilterType BySJout",
                          "--outFilterMultimapNmax 20",
                          "--alignSJDBoverhangMin 1",
                          "--outFilterMismatchNMax 999",
                          "--alignIntronMin 20",
                          "--alignIntronMax 1000000",
                          "--alignMatesGapMax 1000000"]
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
    output_sam = "{}Aligned.out.sam".format(output)
    logfile = "{}.star.log".format(name)
    sort_logfile = "{}.sortconvert.log".format(name)
    temp = "{}.bwa.sort.temp".format(name)
    output_bam = "{}.star.bam".format(name)

    command = ["{}".format(config['star']['bin']),
               "--genomeDir {}".format(config['star']['index']),
               "--runThreadN {}".format(config['star']['num_cores']),
               "--readFilesIn {} {}".format(samples[name]['fastq1'], samples[name]['fastq2']),
               "--outFileNamePrefix {}".format(output),
               "--outReadsUnmapped Fastx"
               ]

    command = add_additional_options(command, config, flags)

    view_cmd = ["{}".format(config['samtools']['bin']),
                "view",
                "-T {}".format(config['reference']),
                "-u",
                "{}".format(output_sam)
                ]

    sort_cmd = ["{}".format(config['samtools']['bin']),
                "sort",
                "-@",
                "{}".format(config['star']['num_cores']),
                "-O",
                "bam",
                "-o",
                "{}".format(output_bam),
                "-T",
                "{}".format(temp),
                "-"
                ]

    sort_command = "{} | {}".format(" ".join(view_cmd), " ".join(sort_cmd))

    job.fileStore.logToMaster("STAR Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    job.fileStore.logToMaster("Sort and Convert Command: {}\n".format(sort_command))
    pipeline.run_and_log_command(sort_command, sort_logfile)

    samples[name]['bam'] = output_bam

    return output_bam


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

    command = ["{}".format(config['star']['bin']),
               "--genomeDir {}".format(config['star']['index']),
               "--runThreadN {}".format(config['star']['num_cores']),
               "--readFilesIn {}".format(samples[name]['fastq1']),
               "--outFileNamePrefix {}".format(output),
               "--outReadsUnmapped Fastx"
               ]

    command = add_additional_options(command, config, flags)

    job.fileStore.logToMaster("STAR Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    aligned_sam = "{}Aligned.out.sam".format(output)
    samples[name]['bam'] = aligned_sam

    return output
