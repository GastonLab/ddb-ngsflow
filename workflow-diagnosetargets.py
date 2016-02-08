#!/usr/bin/env python

# Standard packages
import sys
import argparse

# Third-party packages
from toil.job import Job

# Package methods
from ddb import configuration
from ngsflow import gatk
from ngsflow.utils import utilities


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    # args.logLevel = "INFO"

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    samples = configuration.configure_samples(args.samples_file, config)

    root_job = Job.wrapJobFn(utilities.spawn_batch_jobs)

    for sample in samples:
        diagnose_targets_job = Job.wrapJobFn(gatk.diagnosetargets, config, sample, samples, samples[sample]['bam'],
                                             cores=int(config['gatk']['num_cores']),
                                             memory="{}G".format(config['gatk']['max_mem']))
        root_job.addChild(diagnose_targets_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
