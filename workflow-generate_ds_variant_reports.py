#!/usr/bin/env python

# Standard packages
import sys
import argparse

# Third-party packages
from toil.job import Job

# Package methods
from ngsflow.variation import variation
from ngsflow.utils import configuration
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
    samples = configuration.configure_samples(args.samples_file)

    root_job = Job.wrapJobFn(utilities.spawn_batch_jobs)
    return_files = list()

    sys.stdout.write("Processing samples:\n")
    for sample in samples:
        sys.stdout.write("{}\n".format(sample))

    for sample in samples:
        report_job = Job.wrapJobFn(variation.generate_ds_variant_report, config, sample, samples[sample]['db1'],
                                   samples[sample]['db2'],
                                   cores=1, memory="2G")
        root_job.addChild(report_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
