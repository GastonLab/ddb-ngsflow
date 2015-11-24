__author__ = 'dgaston'

__author__ = 'dgaston'

# Standard packages
import sys
import argparse
import multiprocessing

# Third-party packages
from toil.job import Job

# Package methods
from ngsflow import gatk
from ngsflow import annotation
from ngsflow import read_sample_sheet
from ngsflow.align import bwa
from ngsflow.utils import configuration
from ngsflow.utils import utilities
from ngsflow.variation import variation
from ngsflow.variation import freebayes
from ngsflow.variation import mutect
from ngsflow.variation import platypus
from ngsflow.variation import vardict
from ngsflow.variation import scalpel
from ngsflow.variation import indelminer


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="Input configuration file for samples")
    parser.add_argument('-c', '--configuration', help="Configuration file for various settings")
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()
    # args.logLevel = "INFO"

    sys.stdout.write("Parsing configuration data\n")
    config = configuration.configure_runtime(args.configuration)

    sys.stdout.write("Parsing sample data\n")
    samples = read_sample_sheet.read(args.samples_file)