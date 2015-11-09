__author__ = 'dgaston'

import sys
import subprocess as sub
import multiprocessing

from toil.job import Job


def run_bwa_mem(job, config, sample_id, fastq1, fastq2):
    """Run BWA MEM  and pipe to samtoools to sort and convert to BAM format"""

    # task_desc = "BWA-MEM: %s" % sample_config[sample_id]['name']

    job.fileStore.logToMaster("Running BWA for sample {}\n".format(sample_id))

    output_bam = "{}.bwa.sorted.bam".format(sample_id)
    temp = "{}.bwa.sort.temp".format(sample_id)

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

    logfile = "{}.alignment.log".format(output_bam)
    with open(logfile, "wb") as err:
        sys.stdout.write("Executing {} and writing to logfile %s\n".format(command))
        err.write("Command: {}\n".format(command))
        job.fileStore.logToMaster("BWA Command: {}\n".format(command))
        # p = sub.Popen(command, stdout=sub.PIPE, stderr=err, shell=True)
        # output = p.communicate()
        # code = p.returncode
        # if code:
        #     sys.stdout.write("An error occurred. Please check %s for details\n" % logfile)
        #     sys.stdout.write("%s\n" % output)
        #     sys.stderr.write("An error occurred. Please check %s for details\n" % logfile)
