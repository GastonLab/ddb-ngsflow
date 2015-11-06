__author__ = 'dgaston'

import sys
import argparse
import multiprocessing
import subprocess as sub


def run_bwa_mem(config, fastq1, fastq2, output_bam):
    """Run BWA MEM  and pipe to samtoools to sort and convert to BAM format"""

    # task_desc = "BWA-MEM: %s" % sample_config[sample_id]['name']

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

    # command = "{} | {} | {}".format(" ".join(bwa_cmd), " ".join(view_cmd), " ".join(sort_cmd))
    command = ("""%s mem -t %s -M -v 2 %s %s %s | %s view -b -S -u - | %s sort -@ %s - %s""" %
               (config['bwa'], multiprocessing.cpu_count(), config['reference_genome'], fastq1, fastq2,
                config['samtools'], config['samtools'], multiprocessing.cpu_count(), output_bam))
    logfile = "{}.alignment.log".format(output_bam)
    with open(logfile, "wb") as err:
        sys.stdout.write("Executing {} and writing to logfile {}\n".format(command, logfile))
        err.write("Command: {}\n".format(command))
        p = sub.Popen(command, stdout=sub.PIPE, stderr=err, shell=True)
        output = p.communicate()
        code = p.returncode
        if code:
            sys.stdout.write("An error occurred. Please check %s for details\n" % logfile)
            sys.stdout.write("%s\n" % output)
            sys.stderr.write("An error occurred. Please check %s for details\n" % logfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-1', '--fastq1', help="FastQ1")
    parser.add_argument('-2', '--fastq2', help="FastQ1")
    parser.add_argument('-o', '--output', help="Output BAM file name")

    args = parser.parse_args()

    configuration = dict()
    configuration['reference_genome'] = "/data/Resources/Genomes/Human/GATK-Bundle/2.8/b37/human_g1k_v37.fasta"
    configuration['samtools'] = "samtools"
    configuration['bwa'] = "bwa"

    run_bwa_mem(configuration, args.fastq1, args.fastq2, args.output)