#!/usr/bin/env python

# Standard packages
import sys
import argparse

# Third-party packages
from toil.job import Job

# Package methods
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

    sys.stdout.write("Processing samples:\n")
    for sample in samples:
        sys.stdout.write("{}\n".format(sample))

    sample_results = dict()
    for sample in samples:
        scan_job = Job.wrapJobFn(utilities.read_coverage, config, sample, samples[sample]['vcf'],
                                 cores=1, memory="2G")
        sample_results[sample] = scan_job.rv()
        root_job.addChild(scan_job)

    summarize_job = Job.wrapJobFn(utilities.generate_coverage_summary, config, sample_results)
    root_job.addFollowOn(summarize_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
