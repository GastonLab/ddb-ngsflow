#!/usr/bin/env python

# Standard packages
import sys
import argparse

# Third-party packages
from toil.job import Job

# Package methods
from ddb import configuration
from ngsflow import gatk
from ngsflow import annotation
from ngsflow.align import bwa
from ngsflow.utils import utilities
from ngsflow.variation import variation
from ngsflow.variation import freebayes
from ngsflow.variation import mutect
from ngsflow.variation import platypus
from ngsflow.variation import vardict
from ngsflow.variation import scalpel
from ngsflow.variation import indelminer
from ngsflow.variation import pisces


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    args.logLevel = "INFO"

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    samples = configuration.configure_samples(args.samples_file, config)

    # Workflow Graph definition. The following workflow definition should create a valid Directed Acyclic Graph (DAG)
    root_job = Job.wrapJobFn(utilities.spawn_batch_jobs, cores=1)
    # root_job.addChildJobFn(utilities.run_fastqc, config, samples,
    #                        cores=1,
    #                        memory="{}G".format(config['fastqc']['max_mem']))

    # Per sample jobs
    for sample in samples:
        # Alignment and Refinement Stages
        align_job = Job.wrapJobFn(bwa.run_bwa_mem, config, sample, samples,
                                  cores=int(config['bwa']['num_cores']),
                                  memory="{}G".format(config['bwa']['max_mem']))

        add_job = Job.wrapJobFn(gatk.add_or_replace_readgroups, config, sample, align_job.rv(),
                                cores=1,
                                memory="{}G".format(config['gatk']['max_mem']))

        creator_job = Job.wrapJobFn(gatk.realign_target_creator, config, sample, add_job.rv(),
                                    cores=int(config['gatk']['num_cores']),
                                    memory="{}G".format(config['gatk']['max_mem']))

        realign_job = Job.wrapJobFn(gatk.realign_indels, config, sample, add_job.rv(), creator_job.rv(),
                                    cores=1,
                                    memory="{}G".format(config['gatk']['max_mem']))

        recal_job = Job.wrapJobFn(gatk.recalibrator, config, sample, realign_job.rv(),
                                  cores=int(config['gatk']['num_cores']),
                                  memory="{}G".format(config['gatk']['max_mem']))

        # Create workflow from created jobs
        root_job.addChild(align_job)
        align_job.addChild(add_job)
        add_job.addChild(creator_job)
        creator_job.addChild(realign_job)
        realign_job.addChild(recal_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
