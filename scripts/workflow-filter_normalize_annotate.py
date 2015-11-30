__author__ = 'dgaston'

# Standard packages
import sys
import argparse
import multiprocessing

# Third-party packages
from toil.job import Job

# Package methods
from ngsflow import gatk
from ngsflow import annotation
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

    num_cores = multiprocessing.cpu_count()

    root_job = Job.wrapJobFn(utilities.spawn_batch_jobs)

    for sample in samples:
        on_target_job = Job.wrapJobFn(utilities.bcftools_filter_variants_regions, config, sample,
                                      samples[sample]['vcf'], cores=num_cores, memory="1G")
        gatk_annotate_job = Job.wrapJobFn(gatk.annotate_vcf, config, sample, on_target_job.rv(), samples[sample]['bam'],
                                          cores=num_cores, memory="4G")
        gatk_filter_job = Job.wrapJobFn(gatk.filter_variants, config, sample, gatk_annotate_job.rv(),
                                        cores=1, memory="2G")
        normalization_job = Job.wrapJobFn(utilities.vt_normalization, config, sample, gatk_filter_job.rv(),
                                          cores=1, memory="2G")
        snpeff_job = Job.wrapJobFn(annotation.snpeff, config, sample, normalization_job.rv(),
                                   cores=num_cores, memory="4G")
        gemini_job = Job.wrapJobFn(annotation.gemini, config, sample, snpeff_job.rv(),
                                   cores=num_cores, memory="4G")

        root_job.addChild(on_target_job)
        on_target_job.addChild(gatk_annotate_job)
        gatk_annotate_job.addChild(gatk_filter_job)
        gatk_filter_job.addChild(normalization_job)
        normalization_job.addChild(snpeff_job)
        snpeff_job.addChild(gemini_job)

    Job.Runner.startToil(root_job, args)
