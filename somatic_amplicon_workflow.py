__author__ = 'dgaston'

# Standard packages
import sys
import argparse
import multiprocessing

# Third-party packages
from toil.job import Job

# Package methods
from ngsflow import gatk
from ngsflow import read_sample_sheet
from ngsflow.align import bwa
from ngsflow.utils import utilities
from ngsflow.variation import variation


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

    # Per sample jobs
    for sample in samples:
        # Alignment and Refinement Stages
        align_job = Job.wrapJobFn(bwa.run_bwa_mem, config, sample, samples[sample]['fastq1'], samples[sample]['fastq2'])
        add_job = Job.wrapJobFn(gatk.add_or_replace_readgroups, config, sample, align_job.rv())
        realign_job = Job.wrapJobFn(gatk.realign_indels, config, sample, add_job.rv())
        recal_job = Job.wrapJobFn(gatk.recalibrator, config, sample, realign_job.rv())

        # Variant calling
        spawn_variant_job = Job.wrapJobFn(utilities.spawn_variant_jobs)
        freebayes_job = Job.wrapJobFn()
        mutect_job = Job.wrapJobFn()
        vardict_job = Job.wrapJobFn()
        scalpel_job = Job.wrapJobFn()
        pindel_job = Job.wrapJobFn()
        indelminer_job = Job.wrapJobFn()
        platypus_job = Job.wrapJobFn()

        # Merge results and annotate
        merge_job = Job.wrapJobFn(variation.merge_variant_calls, config, sample, (freebayes_job.rv(), mutect_job.rv(),
                                  vardict_job.rv(), scalpel_job.rv(), pindel_job.rv(), indelminer_job.rv(),
                                  platypus_job.rv()))
        gatk_anno_filter_job = Job.wrapJobFn()
        normalization_job = Job.wrapJobFn()
        snpeff_job = Job.wrapJobFn()
        gemini_job = Job.wrapJobFn()

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
        spawn_variant_job.addChild(pindel_job)
        spawn_variant_job.addChild(indelminer_job)
        spawn_variant_job.addChild(platypus_job)

        spawn_variant_job = spawn_variant_job.encapsulate()
        spawn_variant_job.addChild(merge_job)
        merge_job.addChild(gatk_anno_filter_job)
        gatk_anno_filter_job.addChild(normalization_job)
        normalization_job.addChild(snpeff_job)
        snpeff_job.addChild(gemini_job)

    # Jobs to be executed for a cohort if necessary
    root_job = root_job.encapsulate()

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
