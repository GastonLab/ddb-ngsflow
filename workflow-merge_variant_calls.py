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
from ngsflow import pipeline
from ngsflow.variation import variation


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
    root_job = Job.wrapJobFn(pipeline.spawn_batch_jobs, cores=1)
    # root_job.addChildJobFn(utilities.run_fastqc, config, samples,
    #                        cores=1,
    #                        memory="{}G".format(config['fastqc']['max_mem']))

    # Per sample jobs
    for sample in samples:
        # Need to filter for on target only results somewhere as well
        spawn_variant_job = Job.wrapJobFn(pipeline.spawn_variant_jobs)
        normalization_job1 = Job.wrapJobFn(variation.vt_normalization, config, sample, "mutect",
                                           samples[sample]['mutect'],
                                           cores=1,
                                           memory="{}G".format(config['gatk']['max_mem']))

        normalization_job2 = Job.wrapJobFn(variation.vt_normalization, config, sample, "scalpel",
                                           samples[sample]['scalpel'],
                                           cores=1,
                                           memory="{}G".format(config['gatk']['max_mem']))

        normalization_job3 = Job.wrapJobFn(variation.vt_normalization, config, sample, "freebayes",
                                           samples[sample]['freebayes'],
                                           cores=1,
                                           memory="{}G".format(config['gatk']['max_mem']))

        normalization_job4 = Job.wrapJobFn(variation.vt_normalization, config, sample, "vardict",
                                           samples[sample]['vardict'],
                                           cores=1,
                                           memory="{}G".format(config['gatk']['max_mem']))

        callers = "mutect,scalpel,freebayes,vardict"

        merge_job = Job.wrapJobFn(variation.merge_variant_calls, config, sample, callers, (normalization_job1.rv(),
                                                                                           normalization_job2.rv(),
                                                                                           normalization_job3.rv(),
                                                                                           normalization_job4.rv()))

        # Removed temporarily until config generation script more easily adds in appropriate region files
        # on_target_job = Job.wrapJobFn(utilities.bcftools_filter_variants_regions, config, sample, samples,
        #                               merge_job.rv())

        gatk_annotate_job = Job.wrapJobFn(gatk.annotate_vcf, config, sample, merge_job.rv(), samples[sample]['bam'],
                                          cores=int(config['gatk']['num_cores']),
                                          memory="{}G".format(config['gatk']['max_mem']))

        gatk_filter_job = Job.wrapJobFn(gatk.filter_variants, config, sample, gatk_annotate_job.rv(),
                                        cores=1,
                                        memory="{}G".format(config['gatk']['max_mem']))

        snpeff_job = Job.wrapJobFn(annotation.snpeff, config, sample, gatk_filter_job.rv(),
                                   cores=int(config['snpeff']['num_cores']),
                                   memory="{}G".format(config['snpeff']['max_mem']))

        gemini_job = Job.wrapJobFn(annotation.gemini, config, sample, snpeff_job.rv(),
                                   cores=int(config['gatk']['num_cores']),
                                   memory="{}G".format(config['gemini']['max_mem']))

        # Create workflow from created jobs
        root_job.addChild(spawn_variant_job)

        spawn_variant_job.addChild(normalization_job1)
        spawn_variant_job.addChild(normalization_job2)
        spawn_variant_job.addChild(normalization_job3)
        spawn_variant_job.addChild(normalization_job4)
        spawn_variant_job.addFollowOn(merge_job)

        merge_job.addChild(gatk_annotate_job)
        # on_target_job.addChild(gatk_annotate_job)
        gatk_annotate_job.addChild(gatk_filter_job)
        gatk_filter_job.addChild(snpeff_job)
        snpeff_job.addChild(gemini_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
