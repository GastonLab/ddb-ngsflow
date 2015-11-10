__author__ = 'dgaston'

import os
import time
import multiprocessing
import subprocess as sub


def touch(path):
    with open(path, 'a'):
        os.utime(path, None)


def spawn_batch_jobs(job):
    """
    This is simply a placeholder root job for the workflow
    """

    job.fileStore.logToMaster("Initializing workflow\n")
    touch("Initialization.bam")
    time.sleep(2)


def spawn_variant_jobs(job):
    """
    This is simply a placeholder job to create a node in the graph for spawning
    off the multiple variant callers
    """

    job.fileStore.logToMaster("Spawning all variant calling methods\n")
    touch("Variant_Calling.bam")
    time.sleep(2)


def run_fastqc(job, config, samples):
    """Run FastQC on provided FastQ files"""

    job.fileStore.logToMaster("Running FastQC for all samples\n")

    fastq_files_list = list()
    for sample in samples:
        fastq_files_list.append(samples[sample]['fastq1'])
        fastq_files_list.append(samples[sample]['fastq2'])

    if multiprocessing.cpu_count() <= len(samples):
        num_cores = multiprocessing.cpu_count()
    else:
        num_cores = len(samples)
    fastq_files_string = " ".join(fastq_files_list)
    command = ("{}".format(config['fastqc']),
               "{}".format(fastq_files_string),
               "--extract",
               "-t",
               "{}".format(num_cores))

    job.fileStore.logToMaster("FastQC Command: {}\n".format(command))
    touch("FastQC.bam")
    time.sleep(2)

    # p = sub.Popen(command, stdout=sub.PIPE, stderr=err, shell=True)
    # output = p.communicate()
    # code = p.returncode
    # if code:
    #     sys.stdout.write("An error occurred. Please check %s for details\n" % logfile)
    #     sys.stdout.write("%s\n" % output)
    #     sys.stderr.write("An error occurred. Please check %s for details\n" % logfile)


def generate_fastqc_summary_report(project_config, sample_config, tool_config, resource_config):
    """Parse FastQC summary reports and generate a run-level summary"""

    with open("%s_fastqc_summary.txt" % project_config['project_name'], 'w') as summary_file:
        for sample in sample_config:
            sample_fastq_dirs = list()
            if sample['fastq1']:
                temp = sample['fastq1'].split(".")
                sample_fastq_dirs.append(temp[0])
            if sample['fastq2']:
                temp = sample['fastq2'].split(".")
                sample_fastq_dirs.append(temp[0])
            for dirbase in sample_fastq_dirs:
                with open("./%s_fastqc/summary.txt" % dirbase, "rU") as fastqc_file:
                    for line in fastqc_file.read():
                        summary_file.write(line)


def run_vt_normalization(project_config, sample_config, tool_config, resource_config):
    """Decompose and left normalize variants"""

    instructions = list()

    if project_config['mode'] == "per_sample":
        for sample in sample_config:
            sample['normalized_vcf'] = "%s.normalized.vcf" % sample['name']
            logfile = "%s.vt_norm.log" % sample['name']

            command = ("zless %s | sed 's/ID=AD.Number=./ID=AD,Number=R/' | vt decompose -s - | "
                       "vt normalize -r %s - > %s" %
                       (sample['working_vcf'], resource_config['reference_genome'], sample['normalized_vcf']))

            instructions.append((command, logfile))
            sample['working_vcf'] = sample['normalized_vcf']
    elif project_config['mode'] == "per_cohort":
        project_config['normalized_vcf'] = "%s.normalized.vcf" % project_config['project_name']
        logfile = "%s.vt_norm.log" % project_config['project_name']

        sys.stdout.write("Performing normalization\n")
        command = ("zless %s | sed 's/ID=AD.Number=./ID=AD,Number=R/' | vt decompose -s - | vt normalize -r %s - > %s" %
                   (project_config['working_vcf'], resource_config['reference_genome'],
                    project_config['normalized_vcf']))
        instructions.append((command, logfile))
        project_config['working_vcf'] = project_config['normalized_vcf']
    else:
        sys.stderr.write("ERROR: Mode: %s not supported\n" % project_config['mode'])
        sys.exit()

    sys.stdout.write("Running normalization\n")
    pipe.execute_multiprocess(instructions, int(tool_config['num_cores']))
    sys.stdout.write("Finished normalization\n")


def generate_coverage_report(project_config, sample_config, tool_config, resource_config):
    """Take DiagnoseTargets data and generate a coverage report"""

    project_config['coverage_problems'] = "%s.coverage_issues_report.txt" % project_config['project_name']
    samples_coverage = {"Chr": [], "Start": [], "Stop": [], "Target": []}
    first_pass = True

    for sample in sample_config:
        filter_field = "%s_filter" % sample['name']
        depth_field = "%s_depth" % sample['name']
        low_field = "%s_bp_low" % sample['name']
        zero_field = "%s_bp_zero" % sample['name']

        samples_coverage[filter_field] = []
        samples_coverage[depth_field] = []
        samples_coverage[low_field] = []
        samples_coverage[zero_field] = []

        targeted_regions = pybedtools.BedTool(resource_config['regions'])
        coverage_data = pybedtools.BedTool(sample['diagnose_targets_vcf'])
        intersections = coverage_data.intersect(targeted_regions, loj=True)

        for region in intersections:
            if first_pass:
                samples_coverage['Chr'].append(region.chrom)
                samples_coverage['Start'].append(region.start)
                samples_coverage['Stop'].append(region.stop)
                samples_coverage['Target'].append(region[13])
            reads_data = region[9].split(":")
            samples_coverage[filter_field].append(region[6])
            samples_coverage[depth_field].append(reads_data[0])
            samples_coverage[low_field].append(reads_data[1])
            samples_coverage[zero_field].append(reads_data[2])

        if first_pass:
            first_pass = False

    content = pyexcel.utils.dict_to_array(samples_coverage)
    sheet = pyexcel.Sheet(content)
    sheet.save_as("%s_coverage_results.xlsx" % project_config['project_name'])


def create_output_dirs(project_config, sample_config, tool_config, resource_config):
    """Create appropriate directories for variant calling outputs and intermediate files"""

    os.mkdir("output")
    os.chdir("output")
    for variant_caller in project_config['variant_callers']:
        os.mkdir(variant_caller)

    os.mkdir("alignment")


def bedtools_coverage(project_config, sample_config, tool_config, resource_config):
    """Use bedtools to evaluate on and off target alignments from bam files"""

    instructions = list()

    for sample in sample_config:
        sample['bedtools_coverage'] = "%s.bedtools_coverage.bed" % sample['name']

        logfile = "%s.bedtools_coverage.log" % sample['name']

        command = ("%s coverage -a %s -b %s -hist" % (tool_config['bedtools']['bin'], sample['working_bam'],
                                                      resource_config['regions']))

        instructions.append((command, logfile))

    sys.stdout.write("Running bedtools to evaluate coverage statistics\n")
    pipe.execute_multiprocess(instructions, int(tool_config['bedtools']['num_cores']))
    sys.stdout.write("Finished bedtools\n")


def bcftools_filter_variants_regions(project_config, sample_config, tool_config, resource_config):
    """Use bcftools to filter vcf file to only variants found within the specified regions file"""

    instructions = list()

    if project_config['mode'] is "per_sample":
        for sample in sample_config:
            sample['on_target_vcf'] = "%s.on_target.vcf" % sample['name']
            logfile = "%s.bcftools.filter.log" % sample['name']

            bgzip_and_tabix_vcf(sample['working_vcf'])
            bgzipped_vcf = "%s.gz" % sample['working_vcf']

            command = ("%s isec -T %s %s > %s" % (tool_config['bcftools']['bin'], resource_config['regions'],
                                                  bgzipped_vcf, sample['on_target_vcf']))
            instructions.append((command, logfile))
            sample['working_vcf'] = sample['on_target_vcf']
    elif project_config['mode'] is "per_cohort":
        project_config['on_target_vcf'] = "%s.on_target.vcf" % project_config['project_name']
        logfile = "%s.bcftools.filter.log" % project_config['project_name']

        bgzip_and_tabix_vcf(project_config['working_vcf'])
        bgzipped_vcf = "%s.gz" % project_config['working_vcf']

        command = ("%s isec -T %s %s > %s" % (tool_config['bcftools']['bin'], resource_config['regions'],
                                              bgzipped_vcf, project_config['on_target_vcf']))
        instructions.append((command, logfile))
        project_config['working_vcf'] = project_config['on_target_vcf']
    else:
        sys.stderr.write("ERROR: Mode: %s not supported\n" % project_config['mode'])
        sys.exit()

    sys.stdout.write("Running bcftools isec to filter to only variants in targets\n")
    pipe.execute_multiprocess(instructions, int(tool_config['bcftools']['num_cores']))
    sys.stdout.write("Finished filtering\n")


def bgzip_and_tabix_vcf_instructions(infile):
    """Generate instructions and logfile for bgzip and tabix"""

    bgzip_command = "bgzip -c %s > %s.gz" % (infile, infile)
    bgzip_logfile = "%s.bgzip.log" % infile

    tabix_command = "tabix -p vcf %s.gz" % infile
    tabix_logfile = "%s.tabix.log" % infile

    bgzip_instructions = list()
    bgzip_instructions.append(bgzip_command)
    bgzip_instructions.append(bgzip_logfile)

    tabix_instructions = list()
    tabix_instructions.append(tabix_command)
    tabix_instructions.append(tabix_logfile)

    return bgzip_instructions, tabix_instructions


def bgzip_and_tabix_vcf(infile):
    """Call bgzip and tabix on vcf files"""

    bgzip_instructions, tabix_instructions = bgzip_and_tabix_vcf_instructions(infile)

    code = pipe.run_and_log_command(bgzip_instructions[0], bgzip_instructions[1])
    pipe.check_return_codes(code)

    code2 = pipe.run_and_log_command(tabix_instructions[0], tabix_instructions[1])
    pipe.check_return_codes(code2)
