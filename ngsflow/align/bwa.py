__author__ = 'dgaston'

import sys
import subprocess as sub
import multiprocessing

from toil.job import Job


def run_bwa_mem(job, sample_id, sample_config, tool_config, resource_config):
    """Run BWA MEM  and pipe to samtoools to sort and convert to BAM format"""

    # task_desc = "BWA-MEM: %s" % sample_config[sample_id]['name']

    sample_config[sample_id]['sorted_bam'] = "%s.bwa.sorted.bam" % sample_config[sample_id]['name']

    bwa_cmd = ["{}".format(config['bwa']),
               "mem",
               "-t",
               "{}".format(multiprocessing.cpu_count()),
               "-M",
               "-v",
               "2",
               "{}".format(config['reference_genome']),
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
                "-o",
                "{}".format(output_bam),
                "-O",
                "bam",
                "-"]
    command = "{} | {} | {}".format(" ".join(bwa_cmd), " ".join(view_cmd), " ".join(sort_cmd))

    logfile = "{}.alignment.log".format(output_bam)
    with open(logfile, "wb") as err:
        sys.stdout.write("Executing {} and writing to logfile %s\n".format(command))
        err.write("Command: {}\n".format(command))
        p = sub.Popen(command, stdout=sub.PIPE, stderr=err, shell=True)
        output = p.communicate()
        code = p.returncode
        if code:
            sys.stdout.write("An error occurred. Please check %s for details\n" % logfile)
            sys.stdout.write("%s\n" % output)
            sys.stderr.write("An error occurred. Please check %s for details\n" % logfile)
