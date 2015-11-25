__author__ = 'dgaston'

import sys
import csv
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="Input file with list of input files to process")
    parser.add_argument('-o', '--output', help="Output coverage report")
    parser.add_argument('-t', '--threshold', help="Threshold for counting low coverage", default=500)
    args = parser.parse_args()

    files = list()
    with open(args.input, 'r') as infile:
        reader = csv.reader(infile)
        for row in reader:
            files.append(row[0])

    coverage_stats = dict()
    for infile in files:
        num_low = 0
        num_bases = 0
        with open(infile, 'r') as infh:
            reader = csv.reader(infh, dialect='excel-tab')
            for row in reader:
                num_bases += 1
                if int(row[-1]) < args.threshold:
                    num_low += 1
            percent = (float(num_low) / float(num_bases)) * 100
            sys.stdout.write("Sample {sample}\t"
                             "Number of Bases Below Threshold: {num}\t"
                             "Percentage: {perc}\n".format(sample=infile, num=num_low, perc=percent))
