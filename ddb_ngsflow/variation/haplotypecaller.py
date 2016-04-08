"""
.. module:: haplotypecaller
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling MuTect.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>

"""

from ddb_ngsflow import pipeline


def haplotypecaller_single(job, config, sample, samples, input_bam):
    """Generate gVCF files for a sample using the HaplotypeCaller
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

    gvcf = "{}.haplotypecaller.g.vcf".format(sample)
    logfile = "{}.haplotypecaller_gvcf.log".format(sample)

    command = ("{}".format(config['gatk']['bin']),
               "-T",
               "HaplotypeCaller",
               "-R",
               "{}".format(config['reference']),
               "--dbsnp",
               "{}".format(config['dbsnp']),
               "-I",
               "{}".format(input_bam),
               "-L",
               "{}".format(samples[sample]['regions']),
               "--emitRefConfidence GVCF",
               "--variant_index_type LINEAR",
               "--variant_index_parameter 128000",
               "-o",
               "{}".format(gvcf))

    job.fileStore.logToMaster("HaplotypeCaller Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return gvcf


def joint_variant_calling(job, config, sample, samples):
    """Create a cohort VCF file based on joint calling from gVCF files
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

    vcf = "{}.haplotypecaller.vcf".format(config['project'])
    logfile = "{}.haplotypecaller_gvcf.log".format(config['project'])

    gvcfs = list()
    for sample in samples:
        gvcfs.append("--variant {}.haplotypecaller.g.vcf".format(sample))

    gvcf_string = " ".join(gvcfs)

    command = ("{}".format(config['gatk']['bin']),
               "-T",
               "GenotypeGVCFs",
               "-R",
               "{}".format(config['reference']),
               "{}".format(gvcf_string),
               "-nt 24",
               "-o",
               "{}".format(vcf))

    job.fileStore.logToMaster("GenotypeVCFs Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return vcf
