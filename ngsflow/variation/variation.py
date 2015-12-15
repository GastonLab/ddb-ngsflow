"""
.. module:: variation
   :platform: Unix, OSX
   :synopsis: A module for various methods and utilities for dealing with variation and VCF files.

.. moduleauthor:: Daniel Gaston <daniel.gaston@gmail.com>


"""
import csv
from gemini import GeminiQuery

from ngsflow.utils import utilities
from ngsflow import pipeline


def _var_is_rare(variant_data):
    """Determine if the MAF of the variant is < 1% in all populations"""

    if variant_data['in_esp'] != 0 or variant_data['in_1kg'] != 0 or variant_data['in_exac'] != 0:
        if variant_data['max_aaf_all'] > 0.01:
            return False
        else:
            return True
    else:
        return True


def _var_is_in_cosmic(variant_data):
    """Determine if the variant has a COSMIC identifier"""

    if variant_data['cosmic_ids'] is not None:
        return True
    else:
        return False


def _var_is_in_clinvar(variant_data):
    """Determine if there is ClinVar data for the variant"""

    if variant_data['clinvar_sig'] is not None:
        return True
    else:
        return False


def _var_not_benign(variant_data):
    """Determine if a variant in clinvar is completely benign"""

    if "benign" in variant_data['clinvar_sig']:
        if "pathogenic" in variant_data['clinvar_sig']:
            return True
        else:
            return False
    else:
        return True


def _var_is_protein_effecting(variant_data):
    if variant_data['impact_severity'] != "LOW":
        return True
    else:
        return False


def _run_gemini_query_and_filter(db):
    """Use the GeminiQuery API to filter results based on severity and specific annotations
    :param db: GEMINI database.
    :type db: str.
    :returns:  tuple -- The header line for the requested columns and all rows that pass filters.
    """

    query = "SELECT chrom, start, end, ref, alt, vcf_id, rs_ids, cosmic_ids, filter, qual, qual_depth, depth, " \
            "(gts).(*), (gt_depths).(*), (gt_ref_depths).(*), (gt_alt_depths).(*), " \
            "gene, transcript, exon, codon_change, aa_change, biotype, impact, impact_so, impact_severity, aa_length, " \
            "is_lof, is_conserved, pfam_domain, " \
            "in_omim, clinvar_sig, clinvar_disease_name, clinvar_origin, clinvar_causal_allele, clinvar_dbsource, " \
            "clinvar_dbsource_id, clinvar_on_diag_assay, " \
            "rmsk, in_segdup, strand_bias, rms_map_qual, in_hom_run, num_mapq_zero, num_reads_w_dels, grc, " \
            "gms_illumina, in_cse, " \
            "num_alleles, allele_count, haplotype_score, is_somatic, somatic_score, " \
            "aaf_esp_ea, aaf_esp_aa, aaf_esp_aa, aaf_esp_all, aaf_1kg_amr, aaf_1kg_eas, aaf_1kg_sas, aaf_1kg_afr, " \
            "aaf_1kg_eur, aaf_1kg_all, aaf_exac_all, aaf_adj_exac_all, aaf_adj_exac_afr, aaf_adj_exac_amr, " \
            "aaf_adj_exac_eas, aaf_adj_exac_fin, aaf_adj_exac_nfe, aaf_adj_exac_oth, aaf_adj_exac_sas, " \
            "max_aaf_all, in_esp, in_1kg, in_exac, " \
            "info FROM variants"
    gq = GeminiQuery(db)
    gq.run(query)
    header = gq.header
    passing_rows = []
    print header

    # Filter out variants with minor allele frequencies above the threshold but
    # retain any that are above the threshold but in COSMIC or in ClinVar and not listed as benign.
    for variant_data in gq:
        if _var_is_in_cosmic(variant_data):
            passing_rows.append(variant_data)
            continue
        if _var_is_in_clinvar(variant_data):
            if _var_not_benign(variant_data):
                passing_rows.append(variant_data)
                continue
        if _var_is_rare(variant_data):
            if _var_is_protein_effecting(variant_data):
                passing_rows.append(variant_data)

    return header, passing_rows


def generate_variant_report(job, config, sample, database):
    """Call the GEMINI Query API and generate a text variant report from the provided database
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param database: The GEMINI database to query.
    :type database: str.
    :returns:  str -- The output vcf file name.
    """

    filename = "{}.variant_report.txt".format(sample)
    header, variants = _run_gemini_query_and_filter(database)
    with open(filename, 'w') as outfile:
        outfile.write("{}\n".format(header))
        for variant in variants:
            outfile.write("{}\n".format(variant))


def compare_ds_library_variants(db1, db2):
    """Use the GeminiQuery API to filter results based on severity and specific annotations
    :param db1: GEMINI database from one stranded library for a sample.
    :type db1: str.
    :param db2: GEMINI database for other stranded library from a sample.
    :type db2: str.
    :returns:  tuple -- The header line for the requested columns and all rows that pass filters.
    """


def merge_variant_calls(job, config, sample, vcf_files):
    """Run vcf-isec to merge variant calls from multiple variant callers
    :param config: The configuration dictionary.
    :type config: dict.
    :param sample: sample name.
    :type sample: str.
    :param vcf_files: List of input vcf files for merging.
    :type vcf_files: list.
    :returns:  str -- The output vcf file name.
    """

    files = list()
    for vcf in vcf_files:
        utilities.bgzip_and_tabix_vcf(job, vcf)
        files.append("{}.gz".format(vcf))
    vcf_files_string = " ".join(files)

    merged_vcf = "{}.merged.vcf".format(sample)
    logfile = "{}.merged.log".format(sample)

    # Put this back in after run: .format(config['vcftools_isec']['bin'])
    isec_command = ("vcf-isec",
                    "-f",
                    "-n",
                    "+1",
                    "{}".format(vcf_files_string),
                    ">",
                    "{}".format(merged_vcf))

    job.fileStore.logToMaster("Vcftools intersect Command: {}\n".format(isec_command))
    pipeline.run_and_log_command(" ".join(isec_command), logfile)

    return merged_vcf
