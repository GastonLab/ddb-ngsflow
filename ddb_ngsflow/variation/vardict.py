"""
.. module:: vardict
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling VarDictJava.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>


"""

from ddb_ngsflow import pipeline


def vardict_matched():
    """Run VarDict in matched tumor/normal mode"""

    raise NotImplementedError()


def vardict_single(job, config, name, samples, input_bam):
    """Run VarDict on an an unmatched tumour sample and call somatic variants
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    vardict_vcf = "{}.vardict.vcf".format(name)
    logfile = "{}.vardict.log".format(name)

    vardict = ["{}".format(config['vardict']['bin']),
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
               "-B",
               "{}".format(config['vardict']['num_cores']),
               # "-a", the amplicon flag seems to be creating errors
               # "-F 0", Probably don't need this as duplicates aren't marked and ignoring secondary alignment good
               "-f",
               "{}".format(config['min_alt_af']),
               "-N",
               "{}".format(name),
               "-b",
               "{}".format(input_bam),
               "{}".format(samples[name]['regions'])]

    vardict2vcf = ["{}".format(config['vardict2vcf']['bin']),
                   "-E",
                   "-f",
                   "{}".format(config['min_alt_af']),
                   "-N",
                   "{}".format(name)]

    vcfsort = ["{}".format(config['vcftools_sort']['bin']),
               "-c"]

    command = ("{vardict} | {strandbias} | {vardict2vcf} | "
               "{sort} > {vcf}".format(vardict=" ".join(vardict), strandbias=config['vardict_strandbias']['bin'],
                                       vardict2vcf=" ".join(vardict2vcf), sort=" ".join(vcfsort), vcf=vardict_vcf))

    job.fileStore.logToMaster("VarDict Command: {}\n".format(command))
    pipeline.run_and_log_command(command, logfile)

    return vardict_vcf
