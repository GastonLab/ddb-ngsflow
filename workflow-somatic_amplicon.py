#!/usr/bin/env python

# Standard packages
import sys
import argparse

# Third-party packages
from toil.job import Job

# Package methods
from ngsflow import gatk
from ngsflow import annotation
from ngsflow.align import bwa
from ngsflow.utils import configuration
from ngsflow.utils import utilities
from ngsflow.variation import variation
from ngsflow.variation import freebayes
from ngsflow.variation import mutect
from ngsflow.variation import platypus
from ngsflow.variation import vardict
from ngsflow.variation import scalpel
from ngsflow.variation import indelminer


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
    samples = configuration.configure_samples(args.samples_file)

    # Workflow Graph definition. The following workflow definition should create a valid Directed Acyclic Graph (DAG)
    root_job = Job.wrapJobFn(utilities.spawn_batch_jobs, cores=1)
    # root_job.addChildJobFn(utilities.run_fastqc, config, samples,
    #                        cores=1,
    #                        memory="{}G".format(config['fastqc']['max_mem']))

    # Per sample jobs
    for sample in samples:
        # Alignment and Refinement Stages
        align_job = Job.wrapJobFn(bwa.run_bwa_mem, config, sample, samples[sample]['fastq1'], samples[sample]['fastq2'],
                                  cores=int(config['bwa']['num_cores']),
                                  memory="{}G".format(config['bwa']['max_mem']))

        add_job = Job.wrapJobFn(gatk.add_or_replace_readgroups, config, sample, align_job.rv(),
                                cores=1,
                                memory="{}G".format(config['gatk']['max_mem']))

        # add_job_bam = add_job.rv()

        creator_job = Job.wrapJobFn(gatk.realign_target_creator, config, sample, add_job.rv(),
                                    cores=int(config['gatk']['num_cores']),
                                    memory="{}G".format(config['gatk']['max_mem']))

        realign_job = Job.wrapJobFn(gatk.realign_indels, config, sample, add_job.rv(), creator_job.rv(),
                                    cores=1,
                                    memory="{}G".format(config['gatk']['max_mem']))

        recal_job = Job.wrapJobFn(gatk.recalibrator, config, sample, realign_job.rv(),
                                  cores=int(config['gatk']['num_cores']),
                                  memory="{}G".format(config['gatk']['max_mem']))

        # Variant calling
        spawn_variant_job = Job.wrapJobFn(utilities.spawn_variant_jobs)

        freebayes_job = Job.wrapJobFn(freebayes.freebayes_single, config, sample, recal_job.rv(),
                                      cores=1,
                                      memory="{}G".format(config['freebayes']['max_mem']))

        mutect_job = Job.wrapJobFn(mutect.mutect_single, config, sample, recal_job.rv(),
                                   cores=int(config['mutect']['num_cores']),
                                   memory="{}G".format(config['mutect']['max_mem']))

        vardict_job = Job.wrapJobFn(vardict.vardict_single, config, sample, recal_job.rv(),
                                    cores=int(config['vardict']['num_cores']),
                                    memory="{}G".format(config['vardict']['max_mem']))

        scalpel_job = Job.wrapJobFn(scalpel.scalpel_single, config, sample, recal_job.rv(),
                                    cores=int(config['scalpel']['num_cores']),
                                    memory="{}G".format(config['scalpel']['max_mem']))

        # indelminer_job = Job.wrapJobFn(indelminer.indelminer_single, config, sample, recal_job.rv(),
        #                                cores=1,
        #                                memory="{}G".format(config['indelminer']['max_mem']))

        platypus_job = Job.wrapJobFn(platypus.platypus_single, config, sample, recal_job.rv(),
                                     cores=int(config['platypus']['num_cores']),
                                     memory="{}G".format(config['platypus']['max_mem']))

        # Merge results and annotate
        merge_job = Job.wrapJobFn(variation.merge_variant_calls, config, sample, (freebayes_job.rv(), mutect_job.rv(),
                                  vardict_job.rv(), scalpel_job.rv(), platypus_job.rv()),
                                  cores=1)

        gatk_annotate_job = Job.wrapJobFn(gatk.annotate_vcf, config, sample, merge_job.rv(), recal_job.rv(),
                                          cores=int(config['gatk']['num_cores']),
                                          memory="{}G".format(config['gatk']['max_mem']))

        gatk_filter_job = Job.wrapJobFn(gatk.filter_variants, config, sample, gatk_annotate_job.rv(),
                                        cores=1,
                                        memory="{}G".format(config['gatk']['max_mem']))

        normalization_job = Job.wrapJobFn(utilities.vt_normalization, config, sample, gatk_filter_job.rv(),
                                          cores=1,
                                          memory="{}G".format(config['gatk']['max_mem']))

        snpeff_job = Job.wrapJobFn(annotation.snpeff, config, sample, normalization_job.rv(),
                                   cores=int(config['snpeff']['num_cores']),
                                   memory="{}G".format(config['snpeff']['max_mem']))

        gemini_job = Job.wrapJobFn(annotation.gemini, config, sample, snpeff_job.rv(),
                                   cores=int(config['gatk']['num_cores']),
                                   memory="{}G".format(config['gemini']['max_mem']))

        # Create workflow from created jobs
        root_job.addChild(align_job)
        align_job.addChild(add_job)
        add_job.addChild(creator_job)
        creator_job.addChild(realign_job)
        realign_job.addChild(recal_job)

        recal_job.addChild(spawn_variant_job)

        spawn_variant_job.addChild(freebayes_job)
        spawn_variant_job.addChild(mutect_job)
        spawn_variant_job.addChild(vardict_job)
        spawn_variant_job.addChild(scalpel_job)
        # spawn_variant_job.addChild(indelminer_job)
        spawn_variant_job.addChild(platypus_job)

        # Use Follow On for Merge so it will follow after all of the spawned variant calling jobs
        spawn_variant_job.addFollowOn(merge_job)
        merge_job.addChild(gatk_annotate_job)
        gatk_annotate_job.addChild(gatk_filter_job)
        gatk_filter_job.addChild(normalization_job)
        normalization_job.addChild(snpeff_job)
        snpeff_job.addChild(gemini_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
