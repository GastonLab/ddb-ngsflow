__author__ = 'dgaston'

import yaml


def read(infile):
    """Read a sample input sheet to the workflow in yaml format. Each sample is its own section
        with FastQ or BAM format input files specified"""

    samples = dict()
    with open(infile, 'r') as yaml_fh:
        data = yaml.safe_load(yaml_fh)

    for sample in data:
        sample_dict = dict()
        if 'bam' in data[sample].keys():
            sample_dict['bam'] = data[sample]['bam']
        if 'fastq1' in data[sample].keys():
            sample_dict['fastq1'] = data[sample]['fastq1']
        if 'fastq2' in data[sample].keys():
            sample_dict['fastq2'] = data[sample]['fastq2']
        samples[sample] = sample_dict

    return samples
