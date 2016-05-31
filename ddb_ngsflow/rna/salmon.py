"""
.. module:: salmon
   :platform: Unix, OSX
   :synopsis: A module of methods for working with the salmon RNA-Seq program
   into additional formats.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>


"""

from ddb_ngsflow import pipeline


def salmon_paired(job, config, name, samples):
    """Annotate the specified VCF using snpEff
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param samples: The samples info and config dictionary.
    :type samples: dict.
    :returns:  str -- The output vcf file name.
    """

    output_dir = "{}.salmon.output".format(name)
    logfile = "{}.salmon.log".format(name)

    command = ("{} quant".format(config['salmon']['bin']),
               "-i {}".format(config['salmon']['index']),
               "-l {}".format(samples[name]['library_type']),
               "-p {}".format(config['salmon']['num_cores']),
               "--useVBOpt",
               "--numBootstraps {}".format(config['salmon']['num_bootstraps']),
               "--biasCorrect",
               "--useFSPD",
               "-1 {}".format(samples[name]['fastq1']),
               "-2 {}".format(samples[name]['fastq2']),
               "-o {}".format(output_dir)
               )

    job.fileStore.logToMaster("Salmon Command: {}\n".format(command))
    pipeline.run_and_log_command(command, logfile)

    return output_dir


def salmon_unpaired(job, config, name, samples):
    """Annotate the specified VCF using snpEff
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param samples: The samples info and config dictionary.
    :type samples: dict.
    :returns:  str -- The output vcf file name.
    """

    output_dir = "{}.salmon.output".format(name)
    logfile = "{}.salmon.log".format(name)

    command = ("{} quant".format(config['salmon']['bin']),
               "-i {}".format(config['salmon']['index']),
               "-l {}".format(samples[name]['library_type']),
               "-p {}".format(config['salmon']['num_cores']),
               "--useVBOpt",
               "--numBootstraps {}".format(config['salmon']['num_bootstraps']),
               "--biasCorrect",
               "--useFSPD",
               "-r {}".format(samples[name]['fastq1']),
               "-o {}".format(output_dir)
               )

    job.fileStore.logToMaster("Salmon Command: {}\n".format(command))
    pipeline.run_and_log_command(command, logfile)

    return output_dir


def salmonEM_unpaired(job, config, name, samples):
    """Annotate the specified VCF using snpEff
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param samples: The samples info and config dictionary.
    :type samples: dict.
    :returns:  str -- The output vcf file name.
    """

    output_dir = "{}.salmon_quant".format(name)
    logfile = "{}.salmon.log".format(name)

    command = ("{} quant".format(config['salmon']['bin']),
               "-i {}".format(config['salmon']['index']),
               "-l {}".format(samples[name]['library_type']),
               "-p {}".format(config['salmon']['num_cores']),
               "--numBootstraps {}".format(config['salmon']['num_bootstraps']),
               "--biasCorrect",
               "--useFSPD",
               "-r {}".format(samples[name]['fastq1']),
               "-o {}".format(output_dir)
               )

    job.fileStore.logToMaster("Salmon Command: {}\n".format(command))
    pipeline.run_and_log_command(command, logfile)

    return output_dir
