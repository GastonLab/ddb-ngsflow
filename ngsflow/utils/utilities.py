__author__ = 'dgaston'

import sys
import multiprocessing

from toil.job import Job


def spawn_batch_jobs(job):
    """
    This is simply a placeholder root job for the workflow
    """

    job.fileStore.logToMaster("Initializing workflow\n")


def run_fastqc(job, config, samples):
    """Run FastQC on provided FastQ files"""

    job.fileStore.logToMaster("Running FastQC for all samples\n")

    fastq_files_list = list()
    for sample in samples:
        fastq_files_list.append(samples[sample]['fastq1'])
        fastq_files_list.append(samples[sample]['fastq2'])

    if multiprocessing.cpu_count() <= len(samples):
        num_cores = multiprocessing.cpu_count()
    else:
        num_cores = len(samples)
    fastq_files_string = " ".join(fastq_files_list)
    command = ("{}".format(config['fastqc']),
               "{}".format(fastq_files_string),
               "--extract",
               "-t",
               "{}".format(num_cores))

    job.fileStore.logToMaster("FastQC Command: {}\n".format(command))
    # p = sub.Popen(command, stdout=sub.PIPE, stderr=err, shell=True)
    # output = p.communicate()
    # code = p.returncode
    # if code:
    #     sys.stdout.write("An error occurred. Please check %s for details\n" % logfile)
    #     sys.stdout.write("%s\n" % output)
    #     sys.stderr.write("An error occurred. Please check %s for details\n" % logfile)
