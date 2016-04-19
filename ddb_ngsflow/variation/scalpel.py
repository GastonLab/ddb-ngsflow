"""
.. module:: scalpel
   :platform: Unix, OSX
   :synopsis: A wrapper module for calling scalpel.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>


"""

import os

from ddb_ngsflow import pipeline
from toil.job import JobException


def scalpel_single(job, config, name, samples, input_bam):
    """Run Scalpel on an an unmatched tumour sample and call somatic variants
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param input_bam: The input_bam file name to process.
    :type input_bam: str.
    :returns:  str -- The output vcf file name.
    """

    cwd = os.getcwd()
    output_dir = os.path.join(cwd, "{}-scalpel-output".format(name))
    scalpel_vcf = os.path.join(output_dir, "variants.indel.vcf")
    fixed_vcf = "{}.scalpel.vcf".format(name)
    logfile = "{}.scalpel.log".format(name)
    logfile2 = "{}.scalpel_fix.log".format(name)

    scalpel_command = ("{}".format(config['scalpel']['bin']),
                       "--single",
                       "--intarget",
                       # "--covthr",
                       # "3",
                       # "--lowcov",
                       # "1",
                       "--ref",
                       "{}".format(config['reference']),
                       "--bed",
                       "{}".format(samples[name]['regions']),
                       "--format",
                       "vcf",
                       "--numprocs",
                       "{}".format(config['scalpel']['num_cores']),
                       "--bam",
                       "{}".format(input_bam),
                       "--dir",
                       "{}".format(output_dir))

    fix_sample_name_command = ("cat",
                               "{}".format(scalpel_vcf),
                               "|",
                               "sed",
                               "'s/sample/{}/g'".format(name),
                               ">",
                               "{}".format(fixed_vcf))

    job.fileStore.logToMaster("Scalpel Command: {}\n".format(scalpel_command))
    pipeline.run_and_log_command(" ".join(scalpel_command), logfile)

    job.fileStore.logToMaster("Scalpel Fix Command: {}\n".format(fix_sample_name_command))
    pipeline.run_and_log_command(" ".join(fix_sample_name_command), logfile2)

    file_path = os.path.join(cwd, fixed_vcf)
    if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
        return scalpel_vcf
    else:
        job.fileStore.logToMaster("Scalpel ran into a problem and no output was generated for file {}. Check logfile"
                                  "{} for details\n".format(scalpel_vcf, logfile))
        return JobException("Scalpel ran into a problem and no output was generated for file {}. Check logfile"
                            "{} for details\n".format(scalpel_vcf, logfile))
