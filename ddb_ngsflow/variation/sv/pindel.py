"""
.. module:: pindel
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling pindel.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>


"""

from ngsflow import pipeline


def run_pindel(job, config, sample, input_bam):
    """Run pindel on samples"""

    pindel_config = "{}.pindel_config.txt".format(sample)
    output_dir = "{}_pindel".format(sample)
    output_vcf = "{}.pindel.vcf".format(sample)

    logfile = "{}.pindel.log".format(sample)
    vcf_logfile = "{}.pindel2vcf.log".format(sample)

    with open(pindel_config, 'w') as bam_config:
        bam_config.write("%s %s %s\n" % (input_bam, config['insert_size'], sample))

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
