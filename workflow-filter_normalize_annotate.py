#!/usr/bin/env python

import sys
import argparse

from toil.job import Job

from ddb import configuration
from ngsflow import gatk
from ngsflow import annotation
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

    for sample in samples:
        on_target_job = Job.wrapJobFn(utilities.bcftools_filter_variants_regions, config, sample,
                                      samples[sample]['vcf'], cores=1, memory="1G")
        gatk_annotate_job = Job.wrapJobFn(gatk.annotate_vcf, config, sample, on_target_job.rv(), samples[sample]['bam'],
                                          cores=int(config['gatk']['num_cores']),
                                          memory="{}G".format(config['gatk']['max_mem']))
        gatk_filter_job = Job.wrapJobFn(gatk.filter_variants, config, sample, gatk_annotate_job.rv(),
                                        cores=1, memory="{}G".format(config['gatk']['max_mem']))
        normalization_job = Job.wrapJobFn(utilities.vt_normalization, config, sample, gatk_filter_job.rv(),
                                          cores=1, memory="2G")
        snpeff_job = Job.wrapJobFn(annotation.snpeff, config, sample, normalization_job.rv(),
                                   cores=1, memory="{}G".format(config['snpeff']['max_mem']))
        gemini_job = Job.wrapJobFn(annotation.gemini, config, sample, snpeff_job.rv(),
                                   cores=int(config['snpeff']['num_cores']), memory="2G")

        root_job.addChild(on_target_job)
        on_target_job.addChild(gatk_annotate_job)
        gatk_annotate_job.addChild(gatk_filter_job)
        gatk_filter_job.addChild(normalization_job)
        normalization_job.addChild(snpeff_job)
        snpeff_job.addChild(gemini_job)

    Job.Runner.startToil(root_job, args)
