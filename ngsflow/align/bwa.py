__author__ = 'dgaston'

import sys
import time
import subprocess as sub
import multiprocessing

from toil.job import Job

from ngsflow.utils import utilities
from ngsflow import pipeline


def run_bwa_mem(job, config, sample, fastq1, fastq2):
    """Run BWA MEM  and pipe to samtoools to sort and convert to BAM format"""

    # task_desc = "BWA-MEM: %s" % sample_config[sample_id]['name']

    job.fileStore.logToMaster("Running BWA for sample {}\n".format(sample))

    output_bam = "{}.bwa.sorted.bam".format(sample)
    temp = "{}.bwa.sort.temp".format(sample)
    logfile = "{}.bwa-align.log".format(sample)

    bwa_cmd = ["{}".format(config['bwa']),
               "mem",
               "-t",
               "{}".format(multiprocessing.cpu_count()),
               "-M",
               "-v",
               "2",
               "{}".format(config['reference']),
               "{}".format(fastq1),
               "{}".format(fastq2)]

    view_cmd = ["{}".format(config['samtools']),
                "view",
                "-u",
                "-"]

    sort_cmd = ["{}".format(config['samtools']),
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
    utilities.touch("{}".format(output_bam))
    time.sleep(2)
    # pipeline.run_and_log_command(command, logfile)

    return output_bam
