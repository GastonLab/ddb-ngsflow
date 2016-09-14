"""
.. module:: pindel
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling pindel.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>


"""

from ddb_ngsflow import pipeline


def run_pindel(job, config, name, input_bam):
    """Run Pindel caller for InDel Detection
    :param config: The configuration dictionary.
    :type config: dict.
    :param name: sample name.
    :type name: str..
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    pindel_config = "{}.pindel_config.txt".format(name)
    output_dir = "{}_pindel".format(name)
    output_vcf = "{}.pindel.vcf".format(name)

    logfile = "{}.pindel.log".format(name)
    vcf_logfile = "{}.pindel2vcf.log".format(name)

    with open(pindel_config, 'w') as bam_config:
        bam_config.write("%s %s %s\n" % (input_bam, config['insert_size'], name))

    command = ("{}".format(config['pindel']['bin']),
               "-f",
               "{}".format(config['reference']),
               "-c",
               "ALL",
               "-w",
               "{}".format(config['pindel']['window']),
               "-E",
               "{}".format(config['pindel']['sensitivity']),
               "-T",
               "{}".format(config['pindel']['num_cores']),
               "-o",
               "{}".format(output_dir),
               "-i",
               "{}".format(pindel_config))

    pindel2vcf_command = ("{}".format(config['pindel2vcf']['bin']),
                          "-r",
                          "{}".format(config['reference']),
                          "-R",
                          "{}".format(config['snpeff']['reference']),
                          "-d",
                          "{}".format(config['snpeff']['reference']),
                          "-he",
                          "0.01",
                          "-G",
                          "-P",
                          "{}".format(output_dir),
                          "-v",
                          "{}".format(output_vcf))

    job.fileStore.logToMaster("Pindel Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    job.fileStore.logToMaster("Pindel2vcf Command: {}\n".format(pindel2vcf_command))
    pipeline.run_and_log_command(" ".join(pindel2vcf_command), vcf_logfile)

    return output_vcf
