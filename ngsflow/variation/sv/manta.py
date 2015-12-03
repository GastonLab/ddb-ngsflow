import os

from ngsflow import pipeline


def manta(job, config, sample, input_bam):
    """Run MANTA caller for Structural Variant Detection
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The Manta output vcf file name.
    """

    cwd = os.getcwd()
    manta_output_dir = os.path.join(cwd, "{}-manta-output".format(sample))
    manta_results_dir = os.path.join(manta_output_dir, "results/variants")
    manta_vcf = os.path.join(manta_results_dir, "tumorSV.vcf.gz")
    manta_config_log = "{}.manta_config.log".format(sample)
    manta_log = "{}.manta.log".format(sample)

    manta_config_command = ()
    manta_command = ()

    job.fileStore.logToMaster("Manta Configuration Command: {}\n".format(manta_command))
    pipeline.run_and_log_command(" ".join(manta_config_command), manta_config_log)

    job.fileStore.logToMaster("Manta Command: {}\n".format(manta_command))
    pipeline.run_and_log_command(" ".join(manta_command), manta_log)

    return manta_vcf
