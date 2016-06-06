"""
.. module:: utilities
   :platform: Unix, OSX
   :synopsis: A module of methods for various QC erelated functions. This is a general catch all module for methoods
   not better implemented elsewhere.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>

"""

from ddb_ngsflow import pipeline


def run_fastqc(job, config, samples):
    """Run FastQC on provided FastQ files
    :param config: The configuration dictionary.
    :type config: dict.
    :param samples: Samples dictionary
    :type samples: str.
    """

    job.fileStore.logToMaster("Running FastQC for all samples\n")
    logfile = "fastqc.log"

    fastq_files_list = list()
    for sample in samples:
        fastq_files_list.append(samples[sample]['fastq1'])
        fastq_files_list.append(samples[sample]['fastq2'])

    fastq_files_string = " ".join(fastq_files_list)
    command = ["{}".format(config['fastqc']['bin']),
               "{}".format(fastq_files_string),
               "--extract"]

    job.fileStore.logToMaster("FastQC Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)


def generate_fastqc_summary_report(job, config, samples):
    """Parse FastQC summary reports and generate a run-level summary
    :param config: The configuration dictionary.
    :type config: dict.
    :param samples: Samples dictionary
    :type samples: str.
    """

    job.fileStore.logToMaster("Parsing FastQC results to run-level summary file\n")
    with open("{}_fastqc_summary.txt".format(config['run_name']), 'w') as summary_file:
        for sample in samples:
            sample_fastq_dirs = list()
            if sample['fastq1']:
                temp = sample['fastq1'].split(".")
                sample_fastq_dirs.append(temp[0])
            if sample['fastq2']:
                temp = sample['fastq2'].split(".")
                sample_fastq_dirs.append(temp[0])
            for dirbase in sample_fastq_dirs:
                with open("./%s_fastqc/summary.txt" % dirbase, "rU") as fastqc_file:
                    for line in fastqc_file.read():
                        summary_file.write(line)
