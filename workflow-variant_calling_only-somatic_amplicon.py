#!/usr/bin/env python

# Standard packages
import sys
import argparse

# Third-party packages
from toil.job import Job

# Package methods
from ngsflow.utils import configuration
from ngsflow.utils import utilities
from ngsflow.variation import freebayes
from ngsflow.variation import mutect
from ngsflow.variation import platypus
from ngsflow.variation import vardict
from ngsflow.variation import scalpel
from ngsflow.variation import pisces
from ngsflow.variation.sv import scanindel


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

    # Per sample variant calling jobs
    for sample in samples:
        # freebayes_job = Job.wrapJobFn(freebayes.freebayes_single, config, sample, samples[sample]['bam'],
        #                               cores=1,
        #                               memory="{}G".format(config['freebayes']['max_mem']))
        #
        # mutect_job = Job.wrapJobFn(mutect.mutect_single, config, sample, samples[sample]['bam'],
        #                            cores=1,
        #                            memory="{}G".format(config['mutect']['max_mem']))

        mutect2_job = Job.wrapJobFn(mutect.mutect2_single, config, sample, samples[sample]['bam'],
                                    cores=1,
                                    memory="{}G".format(config['mutect']['max_mem']))
        #
        # vardict_job = Job.wrapJobFn(vardict.vardict_single, config, sample, samples, samples[sample]['bam'],
        #                             cores=int(config['vardict']['num_cores']),
        #                             memory="{}G".format(config['vardict']['max_mem']))
        #
        # scalpel_job = Job.wrapJobFn(scalpel.scalpel_single, config, sample, samples, samples[sample]['bam'],
        #                             cores=int(config['scalpel']['num_cores']),
        #                             memory="{}G".format(config['scalpel']['max_mem']))

        # scanindel_job = Job.wrapJobFn(scanindel.scanindel, config, sample, samples, samples[sample]['bam'],
        #                               cores=int(config['scanindel']['num_cores']),
        #                               memory="{}G".format(config['scanindel']['max_mem']))

        # Platypus defined regions not created for all panels. Need to handle stranded pool libraries
        platypus_job = Job.wrapJobFn(platypus.platypus_single, config, sample, samples, samples[sample]['bam'],
                                     cores=int(config['platypus']['num_cores']),
                                     memory="{}G".format(config['platypus']['max_mem']))

        # Create workflow from created jobs
        # root_job.addChild(freebayes_job)
        # root_job.addChild(mutect_job)
        root_job.addChild(mutect2_job)
        # root_job.addChild(vardict_job)
        # root_job.addChild(scalpel_job)
        # root_job.addChild(scanindel_job)
        root_job.addChild(platypus_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
