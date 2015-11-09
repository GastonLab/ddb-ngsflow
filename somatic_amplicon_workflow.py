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
    args.logLevel = "INFO"

    config = {"bwa": "bwa", "samtools": "samtools", "fastqc": "fastqc",
              "reference": "/data/Resources/Genomes/Human/GATK-Bundle/2.8/b37/human_g1k_v37.fasta"}

    sys.stdout.write("Parsing sample data\n")
    samples = read_sample_sheet.read(args.samples_file)

    # Workflow Graph definition. The following workflow definition should create a valid Directed Acyclic Graph (DAG)
    root_job = Job.wrapJobFn(utilities.spawn_batch_jobs)
    root_job.addChildJobFn(utilities.run_fastqc, config, samples)

    for sample in samples:
        align_job = Job.wrapJobFn(bwa.run_bwa_mem, config, sample, samples[sample]['fastq1'], samples[sample]['fastq2'])
        add_job = Job.wrapJobFn(gatk.add_or_replace_readgroups, config, sample)

        root_job.addChild(align_job)
        align_job.addChild(add_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
