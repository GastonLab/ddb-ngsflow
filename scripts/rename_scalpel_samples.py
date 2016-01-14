#!/usr/bin/env python

import os
import sys
import argparse
import subprocess as sub

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dry', help="Do a dry-run. Print all rename statements but do not execute", type=bool)
    args = parser.parse_args()

    root_dir = os.getcwd()
    output_directories = os.listdir(os.getcwd())

    # Check Scalpel
    for directory in output_directories:
        os.chdir(directory)
        sections = directory.partition("-scalpel-")
        sample = sections[0]

        fixed_vcf = os.path.join(root_dir, "{}.scalpel.vcf".format(sample))
        logfile = "{}.scalpel-rename.log".format(sample)

        sys.stdout.write("Fixing Scalpel VCF in directory {}:\n".format(directory))
        fix_sample_name_command = ("cat",
                                   "variants.indel.vcf",
                                   "|",
                                   "sed",
                                   "'s/sample/{}/g'".format(sample),
                                   ">",
                                   "{}".format(fixed_vcf))
        command = " ".join(fix_sample_name_command)

        with open(logfile, "wb") as err:
            sys.stdout.write("Executing %s and writing to logfile %s\n" % (command, logfile))
            err.write("Command: %s\n" % command)
            p = sub.Popen(command, stdout=sub.PIPE, stderr=err, shell=True)
            output = p.communicate()
            code = p.returncode
            if code:
                raise RuntimeError("An error occurred when executing the commandline: {}. "
                                   "Please check the logfile {} for details\n".format(command, logfile))
        os.chdir(root_dir)
