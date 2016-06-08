"""
.. module:: cufflinks
   :platform: Unix, OSX
   :synopsis: A module of methods for working with the cufflinks RNA-Seq programs
   into additional formats.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>


"""

import os
from ddb_ngsflow import pipeline


def cufflinks(job, config, name, samples):
    """Transcriptome assembly with cufflinks
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param input_bam: The input bam file.
    :type input_bam: str.
    :returns:  str -- The output transcriptome from cufflinks.
    """

    outdir = "{}_cufflinks".format(name)
    logfile = "{}.cufflinks.log".format(name)

    working_dir = os.getcwd()
    path = os.path.join(working_dir, outdir)
    os.mkdir(path)
    os.chdir(path)

    command = ["{}".format(config['cufflinks']['bin']),
               "-g {}".format(config['transcript_reference']),
               "-b {}".format(config['reference']),
               "-u",
               "-p {}".format(config['cufflinks']['num_cores']),
               "--library-type {}".format(),
               "{}".format(samples[name]['bam'])]

    job.fileStore.logToMaster("Cufflinks Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    os.chdir(working_dir)

    return path


def cuffmerge(job, config, name, manifest):
    """Merge assembled cufflinks transcriptomes from all samples
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param manifest: The file name containing all of the transcript assemblies to be merged.
    :type manifest: str.
    :returns:  str -- The merged output transcriptome from cufflinks.
    """

    stats_root = "{}_cuffmerge_stats".format(config['run_id'])
    logfile = "{}.cuffmerge.log".format(name)

    command = ["{}".format(config['cuffmerge']['bin']),
               "-g {}".format(config['transcript_reference']),
               "-s {}".format(config['reference']),
               "-p {}".format(config['cuffmerge']['num_cores']),
               "{}".format(manifest)]

    job.fileStore.logToMaster("Cufflinks Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return stats_root


def cuffquant(job, config, name, samples):
    """Run Cuffquant on all samples
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :returns:  str -- The directory name for the cuffquant results.
    """

    outdir = "{}_cuffquant".format(name)
    logfile = "{}.cuffquant.log".format(name)

    command = ["{}".format(config['cuffquant']['bin']),
               "-b {}".format(config['reference']),
               "-p {}".format(config['cuffquant']['num_cores']),
               "-u",
               "{}".format(config['transcript_reference']),
               "{}".format(samples[name]['bam'])
               ]

    job.fileStore.logToMaster("Cuffquant Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return outdir


def cuffnorm(job, config, name):
    """Run Cuffnorm on cuffquant results form all samples
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :returns:  str -- The transcript quantification file name.
    """

    outdir = "{}_cuffquant".format(name)
    logfile = "{}.cuffquant.log".format(name)

    command = ["{}".format(config['cuffquant']['bin']),
               "-b {}".format(config['reference']),
               "-p {}".format(config['cuffquant']['num_cores']),
               "-u",
               "{}".format(config['transcript_reference']),
               "{}".format(samples[name]['bam'])
               ]

    job.fileStore.logToMaster("Cuffquant Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return outdir
