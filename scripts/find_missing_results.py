#!/usr/bin/env python

import os
import sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input list of sample/library names to search for")
    args = parser.parse_args()

    samples = list()
    with open(args.samples_file, 'r') as samples_file:
        samples = samples_file.read().splitlines()
    root_dir = os.getcwd()

    # Check FinalBams
    os.chdir("./FinalBAMs")
    sys.stdout.write("Missing final aligned BAM files:\n")
    for sample in samples:
        if not os.path.isfile("{}.recalibrated.sorted.bam".format(sample)):
            sys.stdout.write("{}\n".format(sample))
    os.chdir(root_dir)

    # Check MuTect
    os.chdir("./MuTect")
    sys.stdout.write("Missing MuTect VCFs:\n")
    for sample in samples:
        if not os.path.isfile("{}.mutect.vcf".format(sample)):
            sys.stdout.write("{}\n".format(sample))
    os.chdir(root_dir)

    # Check VarDict
    os.chdir("./VarDict")
    sys.stdout.write("Missing VarDict VCFs:\n")
    for sample in samples:
        if not os.path.isfile("{}.vardict.vcf".format(sample)):
            sys.stdout.write("{}\n".format(sample))
    os.chdir(root_dir)

    # Check Scalpel
    os.chdir("./Scalpel")
    sys.stdout.write("Missing Scalpel VCFs:\n")
    for sample in samples:
        if not os.path.isdir("{}-scalpel-output".format(sample)):
            sys.stdout.write("No Folder: {}\n".format(sample))
        else:
            os.chdir("{}-scalpel-output/".format(sample))
            if not os.path.isfile("variants.indel.vcf"):
                sys.stdout.write("No Variants file: {}\n".format(sample))
            os.chdir("..")
    os.chdir(root_dir)

    # Check FreeBayes
    os.chdir("./FreeBayes")
    sys.stdout.write("Missing FreeBayes VCFs:\n")
    for sample in samples:
        if not os.path.isfile("{}.freebayes.vcf".format(sample)):
            sys.stdout.write("{}\n".format(sample))
    os.chdir(root_dir)

    # Check Coverage
    os.chdir("./Coverage")
    sys.stdout.write("Missing DiagnoseTargets VCFs:\n")
    for sample in samples:
        if not os.path.isfile("{}.diagnosetargets.vcf".format(sample)):
            sys.stdout.write("{}\n".format(sample))
    os.chdir(root_dir)

    # Check GEMINI Databases
    # os.chdir("./GEMINI")
    # sys.stdout.write("Missing GEMINI databases:\n")
    # for sample in samples:
    #     if not os.path.isfile("{}.snpEff.GRCh37.75.db".format(sample)):
    #         sys.stdout.write("{}\n".format(sample))
    # os.chdir(root_dir)
