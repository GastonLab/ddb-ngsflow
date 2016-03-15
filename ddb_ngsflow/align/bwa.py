"""
.. module:: bwa
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling BWA.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>


"""
from ddb_ngsflow import pipeline


def run_bwa_mem(job, config, sample, samples):
    """Run GATK's DiagnoseTargets against the supplied region

    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param fastq1: Input FastQ File.
    :type fastq1: str.
    :param fastq2: Input FastQ File.
    :type fastq2: str.
    :returns:  str -- Aligned and sorted BAM file name.

    """

    job.fileStore.logToMaster("Running BWA for sample {}\n".format(sample))

    output_bam = "{}.bwa.sorted.bam".format(sample)
    temp = "{}.bwa.sort.temp".format(sample)
    logfile = "{}.bwa-align.log".format(sample)

    bwa_cmd = ["{}".format(config['bwa']['bin']),
               "mem",
               "-t",
               "{}".format(config['bwa']['num_cores']),
               "-M",
               "-v",
               "2",
               "{}".format(config['reference']),
               "{}".format(samples[sample]['fastq1']),
               "{}".format(samples[sample]['fastq2'])]

    view_cmd = ["{}".format(config['samtools']['bin']),
                "view",
                "-u",
                "-"]

    sort_cmd = ["{}".format(config['samtools']['bin']),
                "sort",
                "-@",
                "{}".format(config['bwa']['num_cores']),
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
