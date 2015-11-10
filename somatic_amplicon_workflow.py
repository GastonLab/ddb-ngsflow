__author__ = 'dgaston'

# Standard packages
import sys
import argparse
import multiprocessing

# Third-party packages
from toil.job import Job

# Package methods
from ngsflow import gatk
from ngsflow import read_sample_sheet
from ngsflow.align import bwa
from ngsflow.utils import utilities


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    # args.logLevel = "INFO"

    config = {"bwa": "bwa", "samtools": "samtools", "fastqc": "fastqc", "picard": "picard", "max_mem": "4",
              "gatk": "/usr/local/bin/GenomeAnalysisTK.jar", "indel1": "indel1,vcf", "indel2": "indel2.vcf",
              "reference": "/data/Resources/Genomes/Human/GATK-Bundle/2.8/b37/human_g1k_v37.fasta"}

    sys.stdout.write("Parsing sample data\n")
    samples = read_sample_sheet.read(args.samples_file)

    # Workflow Graph definition. The following workflow definition should create a valid Directed Acyclic Graph (DAG)
    root_job = Job.wrapJobFn(utilities.spawn_batch_jobs)
    root_job.addChildJobFn(utilities.run_fastqc, config, samples)

    # Per sample jobs
    for sample in samples:
        align_job = Job.wrapJobFn(bwa.run_bwa_mem, config, sample, samples[sample]['fastq1'], samples[sample]['fastq2'])
        add_job = Job.wrapJobFn(gatk.add_or_replace_readgroups, config, sample, align_job.rv())
        realign_job = Job.wrapJobFn(gatk.realign_indels, config, sample, add_job.rv())

        add_job.addChild(realign_job)
        align_job.addChild(add_job)
        root_job.addChild(align_job)

    # Jobs to be executed for a cohort
    root_job = root_job.encapsulate()

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
