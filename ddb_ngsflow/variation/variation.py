"""
.. module:: variation
   :platform: Unix, OSX
   :synopsis: A module for various methods and utilities for dealing with variation and VCF files.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>


"""

import os
import sys
import cyvcf2

from ddb import gemini_interface
from gemini import GeminiQuery
from ddb_ngsflow import pipeline
from ddb import vcf_parsing
from cyvcf2 import VCF
from cyvcf2 import Writer


def _run_gemini_query_and_filter(db, genes):
    """Use the GeminiQuery API to filter results based on severity and specific annotations
    :param db: GEMINI database.
    :type db: str.
    :returns:  tuple -- The header line for the requested columns and all rows that pass filters.
    """

    query = "SELECT chrom, start, end, ref, alt, vcf_id, rs_ids, cosmic_ids, filter, qual, qual_depth, depth, " \
            "gene, transcript, exon, codon_change, aa_change, biotype, impact, impact_so, impact_severity, " \
            "aa_length, is_lof, is_conserved, pfam_domain, in_omim, clinvar_sig, clinvar_disease_name, " \
            "clinvar_origin, clinvar_causal_allele, clinvar_dbsource, clinvar_dbsource_id, clinvar_on_diag_assay, " \
            "rmsk, in_segdup, strand_bias, rms_map_qual, in_hom_run, num_mapq_zero, num_reads_w_dels, grc, " \
            "gms_illumina, in_cse, num_alleles, allele_count, haplotype_score, is_somatic, somatic_score, " \
            "aaf_esp_ea, aaf_esp_aa, aaf_esp_aa, aaf_esp_all, aaf_1kg_amr, aaf_1kg_eas, aaf_1kg_sas, aaf_1kg_afr, " \
            "aaf_1kg_eur, aaf_1kg_all, aaf_exac_all, aaf_adj_exac_all, aaf_adj_exac_afr, aaf_adj_exac_amr, " \
            "aaf_adj_exac_eas, aaf_adj_exac_fin, aaf_adj_exac_nfe, aaf_adj_exac_oth, aaf_adj_exac_sas, " \
            "max_aaf_all, in_esp, in_1kg, in_exac FROM variants"
    # "(gts).(*), (gt_depths).(*), (gt_ref_depths).(*), (gt_alt_depths).(*), " \
    gq = GeminiQuery(db)
    gq.run(query)
    header = gq.header
    passing_rows = []
    print header

    # Filter out variants with minor allele frequencies above the threshold but
    # retain any that are above the threshold but in COSMIC or in ClinVar and not listed as benign.
    for variant_data in gq:
        if genes:
            if not gemini_interface.var_in_gene(variant_data, genes):
                continue
        # Right now removing this. Many benign and synonymous variants are in cosmic
        # if _var_is_in_cosmic(variant_data):
        #     passing_rows.append(variant_data)
        #     continue
        if gemini_interface.var_is_in_clinvar(variant_data):
            # Removed is_benign check temporarily. Some variants not annotated with up to date annotations
            passing_rows.append(variant_data)
            continue
        if gemini_interface.var_is_rare(variant_data):
            if gemini_interface.var_is_protein_effecting(variant_data):
                passing_rows.append(variant_data)

    return header, passing_rows


def generate_variant_report(job, config, name, genes, database):
    """Call the GEMINI Query API and generate a text variant report from the provided database
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param database: The GEMINI database to query.
    :type database: str.
    :returns:  None
    """

    filename = "{}.variant_report.txt".format(name)
    header, variants = _run_gemini_query_and_filter(database, genes)
    with open(filename, 'w') as outfile:
        outfile.write("{}\n".format(header))
        for variant in variants:
            outfile.write("{}\n".format(variant))


def vt_normalization(job, config, sample, caller, input_vcf):
    """Decompose and left normalize variants
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param caller: caller name.
    :type caller: str.
    :param input_vcf: The input_vcf file name to process.
    :type input_vcf: str.
    :returns:  str -- The output vcf file name.
    """

    output_vcf = "{}.{}.normalized.vcf".format(sample, caller)
    logfile = "{}.{}.vt_normalization.log".format(sample, caller)

    normalization = ["zless",
                     "{}".format(input_vcf),
                     "|",
                     "sed",
                     "'s/ID=AD,Number=./ID=AD,Number=R/'",
                     "|",
                     "{}".format(config['vt']['bin']),
                     "decompose",
                     "-s",
                     "-",
                     "|",
                     "{}".format(config['vt']['bin']),
                     "normalize",
                     "-r",
                     "{}".format(config['reference']),
                     "-",
                     ">",
                     "{}".format(output_vcf)]

    job.fileStore.logToMaster("VT Command: {}\n".format(normalization))
    pipeline.run_and_log_command(" ".join(normalization), logfile)

    return output_vcf


def bgzip_tabix_vcf(job, config, sample, caller, input_vcf):
    """BGZip and Tabix an input VCF file
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param caller: caller name.
    :type caller: str.
    :returns:  str -- The output vcf file name.
    """

    bgzip_vcf = "{}.gz".format(input_vcf)
    bgzip_cmd = ["bgzip",
                 "-c",
                 "{}".format(input_vcf),
                 ">",
                 "{}".format(bgzip_vcf)]

    tabix_cmd = ["tabix",
                 "-p",
                 "vcf",
                 "{}".format(bgzip_vcf)]

    bgzip_logfile = "{}.{}.bgzip.log".format(sample, caller)
    tabix_logfile = "{}.{}.tabix.log".format(sample, caller)

    job.fileStore.logToMaster("Bgzip Command: {}\n".format(bgzip_cmd))
    pipeline.run_and_log_command(" ".join(bgzip_cmd), bgzip_logfile)

    job.fileStore.logToMaster("Tabix Command: {}\n".format(tabix_cmd))
    pipeline.run_and_log_command(" ".join(tabix_cmd), tabix_logfile)

    return bgzip_vcf


def add_refcontig_info_header(job, config, sample, caller, input_vcf):
    """Use BCFTools to add contig info from reference to VCF
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param caller: caller name.
    :type caller: str.
    :returns:  str -- The output vcf file name.
    """

    output_vcf = "{}.{}.rehead.vcf.gz".format(sample, caller)
    logfile = "{}.{}.bcftools_rehead.log".format(sample, caller)
    command = ["bcftools-1.9dev",
               "reheader",
               "-f",
               "{}.fai".format(config['reference']),
               "-o",
               "{}".format(output_vcf),
               "{}".format(input_vcf)]

    job.fileStore.logToMaster("BCFTools Reheader Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return output_vcf


def PicardUpdateVCFDict(job, config, sample, caller, input_vcf):
    """Use Picard to add contig info from reference to VCF
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param caller: caller name.
    :type caller: str.
    :returns:  str -- The output vcf file name.
    """

    output_vcf = "{}.{}.rehead.vcf".format(sample, caller)
    logfile = "{}.{}.picard_rehead.log".format(sample, caller)

    command = ["{}".format(config['picard']['bin']),
               "UpdateVcfSequenceDictionary",
               "INPUT={}".format(input_vcf),
               "OUTPUT={}".format(output_vcf),
               "SEQUENCE_DICTIONARY={}".format(config['dict'])]

    job.fileStore.logToMaster("Picard Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile)

    return output_vcf


def filter_low_support_variants(job, config, sample, caller, input_vcf):
    """Filter out very low quality calls from VCFs so they are not included in database
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param caller: caller name.
    :type caller: str.
    :returns:  str -- The output vcf file name.
    """

    output_vcf = "{}.{}.low_support_filtered.vcf".format(sample, caller)

    job.fileStore.logToMaster("Filtering VCF {}\n".format(input_vcf))
    parse_functions = {'mutect': vcf_parsing.parse_mutect_vcf_record,
                       'freebayes': vcf_parsing.parse_freebayes_vcf_record,
                       'vardict': vcf_parsing.parse_vardict_vcf_record,
                       'scalpel': vcf_parsing.parse_scalpel_vcf_record,
                       'platypus': vcf_parsing.parse_platypus_vcf_record,
                       'pindel': vcf_parsing.parse_pindel_vcf_record}

    vcf = VCF(input_vcf)
    writer = Writer(output_vcf, vcf)

    for variant in vcf:
        pass_filter = True
        var_info = parse_functions[caller](variant)
        if float(var_info['Alt_Depth']) < 5.0:
            pass_filter = False
        if pass_filter:
            writer.write_record(variant)

    writer.close()
    return output_vcf


def merge_variant_calls(job, config, sample, callers, vcf_files):
    """Merge variant calls from multiple variant callers
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param callers: Comma-separated list of VCF callers to tag the ensemble output. Must be in same order as vcf_files.
    :type sample: str.
    :param vcf_files: List of input vcf files for merging.
    :type vcf_files: list.
    :returns:  str -- The output vcf file name.
    """

    merged_vcf = "{}.merged.vcf.gz".format(sample)
    uncompressed_vcf = "{}.merged.vcf".format(sample)
    sorted_vcf = "{}.merged.sorted.vcf".format(sample)

    logfile1 = "{}.merging.log".format(sample)
    logfile2 = "{}.uncompress-merging.log".format(sample)
    logfile3 = "{}.merged_sort.log".format(sample)

    vcf_files_string = " ".join(vcf_files)

    command = ["{}".format(config['ensemble']['bin']),
               "ensemble",
               "-c",
               "{}".format(config['ensemble']['num_cores']),
               "--numpass",
               "1",
               "--names",
               "{}".format(callers),
               "{}".format(merged_vcf),
               "{}".format(config['reference']),
               "{}".format(vcf_files_string)]

    command2 = ["bgzip",
                "-cd",
                "{}".format(merged_vcf),
                ">",
                "{}".format(uncompressed_vcf)]

    command3 = ["{}".format(config['picard']['bin']),
                "SortVcf",
                "SEQUENCE_DICTIONARY={}".format(config['dict']),
                "OUTPUT={}".format(sorted_vcf),
                "INPUT={}".format(uncompressed_vcf)]

    sys.stderr.write("Running commands: \n")
    sys.stderr.write("bcbio-variation-recall Command: {}\n".format(command))
    sys.stderr.write("Uncompression Command: {}\n".format(command2))
    sys.stderr.write("Sort Command: {}\n".format(command3))

    job.fileStore.logToMaster("bcbio-variation-recall Command: {}\n".format(command))
    pipeline.run_and_log_command(" ".join(command), logfile1)

    job.fileStore.logToMaster("Uncompression Command: {}\n".format(command2))
    pipeline.run_and_log_command(" ".join(command2), logfile2)

    job.fileStore.logToMaster("Sort Command: {}\n".format(command3))
    pipeline.run_and_log_command(" ".join(command3), logfile3)

    # The Index file created by Picard often causes problems with the GATK
    index_file = "{}.idx".format(sorted_vcf)
    os.remove(index_file)

    return sorted_vcf
