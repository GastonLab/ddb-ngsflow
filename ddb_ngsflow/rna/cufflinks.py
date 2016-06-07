"""
.. module:: cufflinks
   :platform: Unix, OSX
   :synopsis: A module of methods for working with the cufflinks RNA-Seq programs
   into additional formats.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>


"""

import os
from ddb_ngsflow import pipeline


def cufflinks(job, config, name, input_bam):
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

    os.mkdir("./{}".format(outdir))
    os.chdir("./{}".format(outdir))

    command = ["{}".format(config['cufflinks']['bin']),
               "-g {}".format(config['transcript_reference']),
               "-b {}".format(config['reference']),
               "-u",
               "-p {}".format(config['cufflinks']['num_cores']),
               "--library-type {}".format(),
               "{}".format(input_bam)]

    job.fileStore.logToMaster("Cufflinks Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return outdir


def cuffmerge(job, config, name, input_transcripts):
    """Merge assembled cufflinks transcriptomes from all samples
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param input_transcripts: list of individual transcript files.
    :type input_transcripts: list.
    :returns:  str -- The merged output transcriptome from cufflinks.
    """
    pass


def cuffquant(job, config, name):
    pass


def cuffnorm(job, config, name):
    pass
