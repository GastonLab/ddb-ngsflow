#!/usr/bin/env python

import os
import sys
import argparse

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
        fixed_vcf = "{}.scalpel.vcf".format(sample)
        sys.stdout.write("Fixing Scalpel VCF in directory {}:\n".format(directory))
        fix_sample_name_command = ("cat",
                                   "variants.indel.vcf",
                                   "|",
                                   "sed",
                                   "'s/sample/{}/g'".format(sample),
                                   ">",
                                   "{}".format(fixed_vcf))
        os.chdir(root_dir)
