"""
.. module:: mutect
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling MuTect.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>

"""

from ddb_ngsflow import pipeline


def mutect_pon():
    """Run MuTect with a synthetic Panel of Normals"""

    raise NotImplementedError()


def mutect_matched():
    """Run MuTect on paired tumor normal data"""

    raise NotImplementedError()


def mutect_single(job, config, name, samples, input_bam):
    """Run MuTect on an an unmatched tumour sample and call somatic variants
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param samples: samples configuration dictionary
    :type samples: dict
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    mutect_vcf = "{}.mutect.vcf".format(name)
    temp_mutect = "{}.tempmutect.vcf".format(name)

    output_stats = "{}.mutectstats.txt".format(name)
    sample_coverage = "{}.mutectcoverage.wig.txt".format(name)

    mutect_logfile = "{}.mutect.log".format(name)
    subset_log = "{}.mutect_subset.log".format(name)

    mutect_command = ("{}".format(config['mutect']['bin']),
                      "-T",
                      "MuTect",
                      "-R",
                      "{}".format(config['reference']),
                      "--dbsnp",
                      "{}".format(config['dbsnp']),
                      "--cosmic",
                      "{}".format(config['cosmic']),
                      "--enable_extended_output",
                      "-I:tumor",
                      "{}".format(input_bam),
                      "--coverage_file",
                      "{}".format(sample_coverage),
                      "-L",
                      "{}".format(samples[name]['regions']),
                      "-isr",
                      "INTERSECTION",
                      "-im",
                      "ALL",
                      "-dt",
                      "NONE",
                      "-o",
                      "{}".format(output_stats),
                      "-vcf",
                      "{}".format(temp_mutect))

    subset_command = ("cat",
                      "{}".format(temp_mutect),
                      "|",
                      "{}".format(config['vcftools_subset']['bin']),
                      "-e",
                      "-c",
                      "{}".format(name),
                      ">",
                      "{}".format(mutect_vcf))

    job.fileStore.logToMaster("MuTect Command: {}\n".format(mutect_command))
    pipeline.run_and_log_command(" ".join(mutect_command), mutect_logfile)

    job.fileStore.logToMaster("Subset Command: {}\n".format(subset_command))
    pipeline.run_and_log_command(" ".join(subset_command), subset_log)

    return mutect_vcf


def mutect2_single(job, config, name, samples, input_bam):
    """Run MuTect on an an unmatched tumour sample and call somatic variants
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    mutect_vcf = "{}.mutect2.vcf".format(name)
    mutect_logfile = "{}.mutect2.log".format(name)

    mutect_command = ("{}".format(config['gatk3.5']['bin']),
                      "-T",
                      "MuTect2",
                      "-R",
                      "{}".format(config['reference']),
                      "--dbsnp",
                      "{}".format(config['dbsnp']),
                      "--cosmic",
                      "{}".format(config['cosmic']),
                      "-drf DuplicateRead",
                      "-ip 100",
                      "-L",
                      "{}".format(samples[name]['regions']),
                      "-nct",
                      "{}".format(config['gatk3.5']['num_cores']),
                      "-I:tumor",
                      "{}".format(input_bam),
                      "-o",
                      "{}".format(mutect_vcf))

    job.fileStore.logToMaster("MuTect2 Command: {}\n".format(mutect_command))
    pipeline.run_and_log_command(" ".join(mutect_command), mutect_logfile)

    # job.fileStore.logToMaster("Subset Command: {}\n".format(subset_command))
    # pipeline.run_and_log_command(" ".join(subset_command), subset_log)

    return mutect_vcf
