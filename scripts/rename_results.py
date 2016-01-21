#!/usr/bin/env python

import os
import sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input list of sample/library names to search for")
    parser.add_argument('-t', '--tag', help="Tag to append to all samples for renaming")
    parser.add_argument('-d', '--dry', help="Do a dry-run. Print all rename statements but do not execute", type=bool)
    args = parser.parse_args()

    samples = list()
    with open(args.samples_file, 'r') as samples_file:
        samples = samples_file.read().splitlines()
    root_dir = os.getcwd()

    # files = os.listdir(os.getcwd())
    # for filename in files:
    #     filename_sections = filename.partition(".")
    #     if filename_sections[0] in samples:
    #         new = "{sample}-{tag}.{end}".format(tag=args.tag, sample=filename_sections[0], end=filename_sections[2])
    #         sys.stdout.write("Renaming file {filename} to {new}\n".format(filename=filename, new=new))
    #         if not args.dry:
    #             sys.stdout.write("Executing\n")
    #             os.rename(filename, new)

    # # Check FastQs
    # os.chdir("./FastQs")
    # sys.stdout.write("Renaming FastQ Files:\n")
    # files = os.listdir(os.getcwd())
    # for filename in files:
    #     filename_list = filename.split("_")
    #     filename_sections = filename.partition(".")
    #     sample_name = "{}_{}".format(filename_list[0], filename_list[1])
    #     if sample_name in samples:
    #         new = "{sample}-{tag}_{read}.{end}".format(tag=args.tag, sample=sample_name, read=filename_list[3],
    #                                                    end=filename_sections[2])
    #         sys.stdout.write("Renaming file {filename} to {new}\n".format(filename=filename, new=new))
    #         if not args.dry:
    #             sys.stdout.write("Executing\n")
    #             os.rename(filename, new)
    # os.chdir(root_dir)
    #
    # Check FinalBams
    os.chdir("./FinalBAMs")
    sys.stdout.write("Renaming final aligned BAM files:\n")
    files = os.listdir(os.getcwd())
    for filename in files:
        filename_sections = filename.partition(".")
        if filename_sections[0] in samples:
            new = "{sample}-{tag}.{end}".format(tag=args.tag, sample=filename_sections[0], end=filename_sections[2])
            sys.stdout.write("Renaming file {filename} to {new}\n".format(filename=filename, new=new))
            if not args.dry:
                sys.stdout.write("Executing\n")
                os.rename(filename, new)
    os.chdir(root_dir)
    #
    # # Check MuTect
    # os.chdir("./MuTect")
    # sys.stdout.write("Renaming MuTect VCFs:\n")
    # files = os.listdir(os.getcwd())
    # for filename in files:
    #     filename_sections = filename.partition(".")
    #     if filename_sections[0] in samples:
    #         new = "{sample}-{tag}.{end}".format(tag=args.tag, sample=filename_sections[0], end=filename_sections[2])
    #         sys.stdout.write("Renaming file {filename} to {new}\n".format(filename=filename, new=new))
    #         if not args.dry:
    #             sys.stdout.write("Executing\n")
    #             os.rename(filename, new)
    # os.chdir(root_dir)
    # #
    #
    # # Check VarDict
    # os.chdir("./VarDict")
    # sys.stdout.write("Renaming VarDict VCFs:\n")
    # files = os.listdir(os.getcwd())
    # for filename in files:
    #     filename_sections = filename.partition(".")
    #     if filename_sections[0] in samples:
    #         new = "{sample}-{tag}.{end}".format(tag=args.tag, sample=filename_sections[0], end=filename_sections[2])
    #         sys.stdout.write("Renaming file {filename} to {new}\n".format(filename=filename, new=new))
    #         if not args.dry:
    #             sys.stdout.write("Executing\n")
    #             os.rename(filename, new)
    # os.chdir(root_dir)
    #
    # Check FreeBayes
    os.chdir("./FreeBayes")
    sys.stdout.write("Renaming FreeBayes VCFs:\n")
    files = os.listdir(os.getcwd())
    for filename in files:
        filename_sections = filename.partition(".")
        if filename_sections[0] in samples:
            new = "{sample}-{tag}.{end}".format(tag=args.tag, sample=filename_sections[0], end=filename_sections[2])
            sys.stdout.write("Renaming file {filename} to {new}\n".format(filename=filename, new=new))
            if not args.dry:
                sys.stdout.write("Executing\n")
                os.rename(filename, new)
    os.chdir(root_dir)
    #
    # # Check GEMINI Databases
    # os.chdir("./GEMINI")
    # sys.stdout.write("Renaming GEMINI databases:\n")
    # files = os.listdir(os.getcwd())
    # for filename in files:
    #     filename_sections = filename.partition(".")
    #     if filename_sections[0] in samples:
    #         new = "{sample}-{tag}.{end}".format(tag=args.tag, sample=filename_sections[0], end=filename_sections[2])
    #         sys.stdout.write("Renaming file {filename} to {new}\n".format(filename=filename, new=new))
    #         if not args.dry:
    #             sys.stdout.write("Executing\n")
    #             os.rename(filename, new)
    # os.chdir(root_dir)
    #
    # # Check Scalpel
    # os.chdir("./Scalpel")
    # sys.stdout.write("Renaming Scalpel VCFs:\n")
    # files = os.listdir(os.getcwd())
    # for filename in files:
    #     filename_sections = filename.partition("-scalpel-")
    #     if filename_sections[0] in samples:
    #         new = "{sample}-{tag}-scalpel-{end}".format(tag=args.tag, sample=filename_sections[0],
    #                                                     end=filename_sections[2])
    #         if os.path.isdir(new):
    #             # sys.stderr.write("Directory {} exists. Skipping\n".format(new))
    #             continue
    #         else:
    #             sys.stdout.write("Renaming file {filename} to {new}\n".format(filename=filename, new=new))
    #             if not args.dry:
    #                 sys.stdout.write("Executing\n")
    #                 os.rename(filename, new)
    # os.chdir(root_dir)
