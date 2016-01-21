#!/usr/bin/env python

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input list of sample/library names to include in the config file")
    parser.add_argument('-o', '--output_file', help="Filename for output config file")
    args = parser.parse_args()

    samples = list()
    with open(args.samples_file, 'r') as samples_file:
        samples = samples_file.read().splitlines()

    with open(args.output_file, 'w') as output_file:
        for sample in samples:
            output_file.write("[{name}]\n".format(name=sample))
            output_file.write("fastq1: ./FastQs/{name}_L001_R1_001.fastq.gz\n".format(name=sample))
            output_file.write("fastq2: ./FastQs/{name}_L001_R2_001.fastq.gz\n".format(name=sample))
            output_file.write("bam: ./FinalBAMs/{name}.recalibrated.sorted.bam\n".format(name=sample))
            output_file.write("mutect: ./MuTect/{name}.mutect.vcf\n".format(name=sample))
            output_file.write("scalpel: ./Scalpel/{name}.scalpel.vcf\n".format(name=sample))
            output_file.write("freebayes: ./FreeBayes/{name}.freebayes.vcf\n".format(name=sample))
            output_file.write("vardict: ./VarDict/{name}.vardict.vcf\n".format(name=sample))
            output_file.write("db: ./GEMINI/{name}.snpEff.GRCh37.75.db\n".format(name=sample))
            output_file.write("\n")
