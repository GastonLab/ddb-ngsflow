"""
.. module:: bwa
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling BWA.

.. moduleauthor:: Daniel Gaston <daniel.gaston@gmail.com>


"""

__author__ = 'dgaston'

import multiprocessing

from ngsflow import pipeline


def run_bwa_mem(job, config, sample, fastq1, fastq2):
    """
    Run BWA MEM  and pipe to samtoools to sort and convert to BAM format

    Return ``bam_file`` upon successful completion of BWA

    """

    job.fileStore.logToMaster("Running BWA for sample {}\n".format(sample))

    output_bam = "{}.bwa.sorted.bam".format(sample)
    temp = "{}.bwa.sort.temp".format(sample)
    logfile = "{}.bwa-align.log".format(sample)

    bwa_cmd = ["{}".format(config['bwa']['bin']),
               "mem",
               "-t",
               "{}".format(multiprocessing.cpu_count()),
               "-M",
               "-v",
               "2",
               "{}".format(config['reference']),
               "{}".format(fastq1),
               "{}".format(fastq2)]

    view_cmd = ["{}".format(config['samtools']['bin']),
                "view",
                "-u",
                "-"]

    sort_cmd = ["{}".format(config['samtools']['bin']),
                "sort",
                "-@",
                "{}".format(multiprocessing.cpu_count()),
                "-O",
                "bam",
                "-o",
                "{}".format(output_bam),
                "-T",
                "{}".format(temp),
                "-"]
    command = "{} | {} | {}".format(" ".join(bwa_cmd), " ".join(view_cmd), " ".join(sort_cmd))

    job.fileStore.logToMaster("BWA Command: {}\n".format(command))
    pipeline.run_and_log_command(command, logfile)

    return output_bam
