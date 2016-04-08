import sys

from ddb_ngsflow import pipeline


def run_delly2_single(job, config, name, input_bam):
    """Run delly2 for structural variant detection. As delly2 is parallelized on the level of samples,
    we use a single-threaded version
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The merged Delly output vcf file name.
    """

    delly_vcfs = list()
    delly_command_core = ("{}".format(config['delly']['bin']),
                          "-x",
                          "{}".format(config['delly']['exclude']),
                          "-g",
                          "{}".format(config['reference']))

    for mut_type in ["DEL", "DUP", "TRA", "INV"]:
        output_vcf = "{sample}.{type}.vcf".format(sample=name, type=mut_type)
        logfile = "{sample}.{type}.log".format(sample=name, type=mut_type)

        delly_vcfs.append(output_vcf)

        delly_command = list()
        delly_command.append(delly_command_core)
        delly_command.append("-t",
                             "{}".format(mut_type),
                             "-o",
                             "{}".format(output_vcf),
                             "{}".format(input_bam))

        job.fileStore.logToMaster("Running Delly: {}\n".format(delly_command))
        pipeline.run_and_log_command(" ".join(delly_command), logfile)

    job.fileStore.logToMaster("Merging delly output with command: {}\n".format(merge_command))
    pipeline.run_and_log_command(" ".join(merge_command), merge_log)

    return merged_vcf
