"""
.. module:: cufflinks
   :platform: Unix, OSX
   :synopsis: A module of methods for working with the cufflinks RNA-Seq programs
   into additional formats.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>


"""

import os
import sys
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
    try:
        os.mkdir(path)
    except:
        sys.stderr.write("Directory {} already exists. Not creating...\n".format(path))

    os.chdir(path)

    samples[name]['cufflinks_assembly'] = os.path.join(path, "transcripts.gtf")

    command = ["{}".format(config['cufflinks']['bin']),
               "-g {}".format(config['transcript_reference']),
               # "-b {}".format(config['reference']),
               "-u",
               "-p {}".format(config['cufflinks']['num_cores']),
               "--library-type {}".format(samples[name]['cufflinks_lib']),
               "{}".format(samples[name]['bam'])]

    if not os.path.isfile("transcripts.gtf"):
        job.fileStore.logToMaster("Cufflinks Command: {}\n".format(command))
        pipeline.run_and_log_command(" ".join(command), logfile)
    else:
        job.fileStore.logToMaster("Cufflinks appears to have already executed for {}. Skipping...\n".format(name))

    os.chdir(working_dir)

    return path


def cuffmerge(job, config, name, samples):
    """Merge assembled cufflinks transcriptomes from all samples
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param samples: Samples config data
    :type samples: dict.
    :returns:  str -- The merged output transcriptome from cufflinks.
    """

    stats_root = "{}_cuffmerge_stats".format(config['run_id'])
    logfile = "{}.cuffmerge.log".format(name)

    manifest_file = "transcript_assemblies.txt"
    with open(manifest_file, 'w') as manifest:
        for sample in samples:
            manifest.write("{}\n".format(samples[sample]['cufflinks_assembly']))

    command = ["{}".format(config['cuffmerge']['bin']),
               "-g {}".format(config['transcript_reference']),
               "-s {}".format(config['reference']),
               "-p {}".format(config['cuffmerge']['num_cores']),
               "{}".format(manifest)]

    job.fileStore.logToMaster("Cufflinks Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    pwd = os.getcwd()
    config['transcript_reference'] = os.path.join(pwd, "merged.gtf")

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
               # "-b {}".format(config['reference']),
               "-p {}".format(config['cuffquant']['num_cores']),
               "-u",
               "{}".format(config['transcript_reference']),
               "{}".format(samples[name]['bam'])
               ]

    job.fileStore.logToMaster("Cuffquant Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return outdir


def cuffnorm(job, config, name, samples):
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
