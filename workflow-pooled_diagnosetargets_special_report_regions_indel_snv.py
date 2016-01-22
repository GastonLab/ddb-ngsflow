#!/usr/bin/env python

# Standard packages
import sys
import argparse

# Third-party packages
from toil.job import Job

# Package methods
from ngsflow import gatk
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
    samples = configuration.configure_samples(args.samples_file, config)

    root_job = Job.wrapJobFn(utilities.spawn_batch_jobs)

    for sample in samples:
        diagnose_snv_targets_job = Job.wrapJobFn(gatk.diagnose_pooled_targets, config, sample, "snv_regions",
                                                 samples, samples[sample]['bam1'], samples[sample]['bam2'],
                                                 cores=int(config['gatk']['num_cores']),
                                                 memory="{}G".format(config['gatk']['max_mem']))

        diagnose_indel_targets_job = Job.wrapJobFn(gatk.diagnose_pooled_targets, config, sample, "indel_regions",
                                                   samples, samples[sample]['bam1'], samples[sample]['bam2'],
                                                   cores=int(config['gatk']['num_cores']),
                                                   memory="{}G".format(config['gatk']['max_mem']))
        root_job.addChild(diagnose_snv_targets_job)
        root_job.addChild(diagnose_indel_targets_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
