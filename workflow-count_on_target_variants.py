__author__ = 'dgaston'

# Standard packages
import sys
import argparse

# Third-party packages
from toil.job import Job

# Package methods
from ngsflow import gatk
from ngsflow.utils import configuration
from ngsflow.utils import utilities


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
    samples = configuration.configure_samples(args.samples_file)

    root_job = Job.wrapJobFn(utilities.spawn_batch_jobs)
    return_files = list()

    sys.stdout.write("Processing samples:\n")
    for sample in samples:
        sys.stdout.write("{}\n".format(sample))

    for sample in samples:
        diagnose_targets_job = Job.wrapJobFn(gatk.diagnosetargets, config, sample, samples[sample]['vcf'],
                                             cores=1, memory="2G")
        return_files.append(diagnose_targets_job.rv())
        root_job.addChild(diagnose_targets_job)

    create_spreadsheet_job = Job.wrapJobFn(utilities.generate_coverage_report, config, return_files)
    root_job.addFollowOn(create_spreadsheet_job)

    # Start workflow execution
    Job.Runner.startToil(root_job, args)
