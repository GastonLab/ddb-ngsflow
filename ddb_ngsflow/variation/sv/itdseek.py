"""
.. module:: freebayes
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling ScanIndel.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>
"""

from ddb_ngsflow import pipeline


def run_flt3_itdseek(job, config, name, samples):
    """Run ITDseek without a matched normal sample
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str.
    :param samples: The samples configuration dictionary.
    :type config: dict.
    :returns:  str -- The output vcf file name.
    """

    itdseek_vcf = "{}.flt3.itdseek.vcf".format(name)
    itdseek_logfile = "{}.flt3.itdseek.log".format(name)

    itdseek_command = ["{}".format(config['itdseek']['bin']),
                       "{}.rg.sorted.bam".format(name),
                       "{}".format(config['reference']),
                       "{}".format(config['samtools-0.19']['bin']),
                       ">",
                       "{}".format(itdseek_vcf)]

    job.fileStore.logToMaster("ITDSeek Command: {}\n".format(itdseek_command))
    pipeline.run_and_log_command(" ".join(itdseek_command), itdseek_logfile)

    return itdseek_vcf
