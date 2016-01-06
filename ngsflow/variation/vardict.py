"""
.. module:: vardict
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling VarDictJava.

.. moduleauthor:: Daniel Gaston <daniel.gaston@gmail.com>


"""

from ngsflow import pipeline


def vardict_matched():
    """Run VarDict in matched tumor/normal mode"""

    raise NotImplementedError()


def vardict_single(job, config, sample, samples, input_bam):
    """Run VarDict on an an unmatched tumour sample and call somatic variants
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    vardict_vcf = "{}.vardict.vcf".format(sample)
    logfile = "{}.vardict.log".format(sample)

    vardict = ("{}".format(config['vardict']['bin']),
               "-G",
               "{}".format(config['reference']),
               "-z",
               "-c",
               "1",
               "-S",
               "2",
               "-E",
               "3",
               "-g",
               "4",
               "-th",
               "{}".format(config['vardict']['num_cores']),
               "-f",
               "{}".format(config['min_alt_af']),
               "-N",
               "{}".format(sample),
               "-b",
               "{}".format(input_bam),
               "{}".format(samples[sample]['regions']))

    vardict2vcf = ("{}".format(config['vardict2vcf']['bin']),
                   "-E",
                   "-f",
                   "{}".format(config['min_alt_af']),
                   "-N",
                   "{}".format(sample))

    command = ("{vardict} | {strandbias} | {vardict2vcf} > {vcf}".format(vardict=" ".join(vardict),
                                                                         strandbias=config['vardict_strandbias']['bin'],
                                                                         vardict2vcf=" ".join(vardict2vcf),
                                                                         vcf=vardict_vcf))

    job.fileStore.logToMaster("VarDict Command: {}\n".format(command))
    pipeline.run_and_log_command(command, logfile)

    return vardict_vcf
