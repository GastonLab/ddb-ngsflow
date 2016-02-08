#!/usr/bin/env python

# Standard packages
import sys
import argparse

# Third-party packages
from toil.job import Job

# Package methods
from ddb import configuration
from ngsflow.variation import variation
from ngsflow import pipeline


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    parser.add_argument('-g', '--genes', help="Comma-separated list of genes to select")
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    # args.logLevel = "INFO"

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    samples = configuration.configure_samples(args.samples_file, config)

    root_job = Job.wrapJobFn(pipeline.spawn_batch_jobs)

    genes = list()
    if args.genes:
        genes = args.genes.split(',')

    for sample in samples:
        report_job = Job.wrapJobFn(variation.generate_variant_report, config, sample, genes, samples[sample]['db'],
                                   cores=1, memory="2G")
        root_job.addChild(report_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
