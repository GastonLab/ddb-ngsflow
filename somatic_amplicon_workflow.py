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
from ngsflow import read_sample_sheet
from ngsflow.align import bwa
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
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    # args.logLevel = "INFO"

    config = {"bwa": "bwa", "samtools": "samtools", "fastqc": "fastqc", "picard": "picard", "max_mem": "4",
              "gatk": "/usr/local/bin/GenomeAnalysisTK.jar", "indel1": "indel1,vcf", "indel2": "indel2.vcf",
              "reference": "/data/Resources/Genomes/Human/GATK-Bundle/2.8/b37/human_g1k_v37.fasta",
              "dbsnp": "/data/Resources/Genomes/Human/GATK-Bundle/2.8/b37/dbsnp_138.b37.vcf"}

    sys.stdout.write("Parsing sample data\n")
    samples = read_sample_sheet.read(args.samples_file)

    # Workflow Graph definition. The following workflow definition should create a valid Directed Acyclic Graph (DAG)
    root_job = Job.wrapJobFn(utilities.spawn_batch_jobs)
    root_job.addChildJobFn(utilities.run_fastqc, config, samples)

    num_cores = multiprocessing.cpu_count()

    # Per sample jobs
    for sample in samples:
        # Alignment and Refinement Stages
        align_job = Job.wrapJobFn(bwa.run_bwa_mem, config, sample, samples[sample]['fastq1'], samples[sample]['fastq2'],
                                  cores=num_cores, memory="4G")
        add_job = Job.wrapJobFn(gatk.add_or_replace_readgroups, config, sample, align_job.rv(),
                                cores=1, memory="4G")
        creator_job = Job.wrapJobFn(gatk.realign_target_creator, config, sample, add_job.rv(),
                                    cores=num_cores, memory="4G")
        realign_job = Job.wrapJobFn(gatk.realign_indels, config, sample, add_job.rv(), creator_job.rv(),
                                    cores=1, memory="4G")
        recal_job = Job.wrapJobFn(gatk.recalibrator, config, sample, realign_job.rv(),
                                  cores=num_cores, memory="4G")

        # Variant calling
        spawn_variant_job = Job.wrapJobFn(utilities.spawn_variant_jobs)

        freebayes_job = Job.wrapJobFn(freebayes.freebayes_single, config, sample, recal_job.rv(),
                                      cores=1, memory="4G")
        mutect_job = Job.wrapJobFn(mutect.mutect_single, config, sample, recal_job.rv(),
                                   cores=num_cores, memory="4G")
        vardict_job = Job.wrapJobFn(vardict.vardict_single, config, sample, recal_job.rv(),
                                    cores=num_cores, memory="5G")
        scalpel_job = Job.wrapJobFn(scalpel.scalpel_single, config, sample, recal_job.rv(),
                                    cores=num_cores, memory="4G")
        indelminer_job = Job.wrapJobFn(indelminer.indelminer_single, config, sample, recal_job.rv(),
                                       cores=1, memory="5G")
        platypus_job = Job.wrapJobFn(platypus.platypus_single, config, sample, recal_job.rv(),
                                     cores=num_cores, memory="4G")

        # Merge results and annotate
        merge_job = Job.wrapJobFn(variation.merge_variant_calls, config, sample, (freebayes_job.rv(), mutect_job.rv(),
                                  vardict_job.rv(), scalpel_job.rv(), indelminer_job.rv(), platypus_job.rv()),
                                  cores=1)
        gatk_annotate_job = Job.wrapJobFn(gatk.annotate_vcf, config, sample, merge_job.rv(), recal_job.rv(),
                                          cores=num_cores, memory="4G")
        gatk_filter_job = Job.wrapJobFn(gatk.filter_variants, config, sample, gatk_annotate_job.rv(),
                                        cores=1, memory="2G")
        normalization_job = Job.wrapJobFn(utilities.vt_normalization, config, sample, gatk_filter_job.rv(),
                                          cores=1, memory="2G")
        snpeff_job = Job.wrapJobFn(annotation.snpeff, config, sample, normalization_job.rv(),
                                   cores=num_cores, memory="4G")
        gemini_job = Job.wrapJobFn(annotation.gemini, config, sample, snpeff_job.rv(),
                                   cores=num_cores, memory="4G")

        # Create workflow from created jobs
        root_job.addChild(align_job)
        align_job.addChild(add_job)
        add_job.addChild(realign_job)
        realign_job.addChild(recal_job)

        recal_job.addChild(spawn_variant_job)
        spawn_variant_job.addChild(freebayes_job)
        spawn_variant_job.addChild(mutect_job)
        spawn_variant_job.addChild(vardict_job)
        spawn_variant_job.addChild(scalpel_job)
        spawn_variant_job.addChild(indelminer_job)
        spawn_variant_job.addChild(platypus_job)

        spawn_variant_job = spawn_variant_job.encapsulate()
        spawn_variant_job.addChild(merge_job)
        merge_job.addChild(gatk_annotate_job)
        gatk_annotate_job.addChild(gatk_filter_job)
        gatk_filter_job.addChild(normalization_job)
        normalization_job.addChild(snpeff_job)
        snpeff_job.addChild(gemini_job)

    # Jobs to be executed for a cohort if necessary
    root_job = root_job.encapsulate()

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
