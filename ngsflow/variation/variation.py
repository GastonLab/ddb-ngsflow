"""
.. module:: variation
   :platform: Unix, OSX
   :synopsis: A module for various methods and utilities for dealing with variation and VCF files.

.. moduleauthor:: Daniel Gaston <daniel.gaston@gmail.com>


"""

__author__ = 'dgaston'

import time

from gemini import GeminiQuery

from ngsflow.utils import utilities
from ngsflow import pipeline


# Can simplify this later with max_aaf_all after updated to 0.18
def var_is_rare(variant_data):
    """Determine if the MAF of the variant is < 1% in all populations"""

    if variant_data['in_esp'] != 0 or variant_data['in_1kg'] != 0 or variant_data['in_exac'] != 0:
        if variant_data['aaf_esp_ea'] > 0.01:
            return False
        elif variant_data['aaf_esp_aa'] > 0.01:
            return False
        elif variant_data['aaf_1kg_amr'] > 0.01:
            return False
        elif variant_data['aaf_1kg_eas'] > 0.01:
            return False
        elif variant_data['aaf_1kg_sas'] > 0.01:
            return False
        elif variant_data['aaf_1kg_afr'] > 0.01:
            return False
        elif variant_data['aaf_1kg_eur'] > 0.01:
            return False
        elif variant_data['aaf_adj_exac_afr'] > 0.01:
            return False
        elif variant_data['aaf_adj_exac_amr'] > 0.01:
            return False
        elif variant_data['aaf_adj_exac_eas'] > 0.01:
            return False
        elif variant_data['aaf_adj_exac_fin'] > 0.01:
            return False
        elif variant_data['aaf_adj_exac_nfe'] > 0.01:
            return False
        elif variant_data['aaf_adj_exac_oth'] > 0.01:
            return False
        elif variant_data['aaf_adj_exac_sas'] > 0.01:
            return False
        else:
            return True
    else:
        return True


def var_is_in_cosmic(variant_data):
    """Determine if the variant has a COSMIC identifier"""

    if variant_data['cosmic_ids'] is not None:
        return True
    else:
        return False


def var_is_in_clinvar(variant_data):
    """Determine if there is ClinVar data for the variant"""

    if variant_data['clinvar_sig'] is not None:
        return True
    else:
        return False


def var_is_protein_effecting(variant_data):
    if variant_data['impact_severity'] != "LOW":
        return True
    else:
        return False


# Need to add max_aaf_all later when gemini updated to 0.18
def run_gemini_query_and_filter(db):
    """Fetch all variants from a specified GEMINI database and filter"""

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
            "in_esp, in_1kg, in_exac, " \
            "info FROM variants"
    gq = GeminiQuery(db)
    gq.run(query)
    header = gq.header
    passing_rows = []
    print header

    # Filter out variants with minor allele frequencies above the threshold but
    # retain any that are above the threshold but in COSMIC
    for variant_data in gq:
        if var_is_in_cosmic(variant_data) or var_is_in_clinvar(variant_data):
            passing_rows.append(variant_data)
            continue
        if var_is_rare(variant_data):
            if var_is_protein_effecting(variant_data):
                passing_rows.append(variant_data)

    return header, passing_rows


def merge_variant_calls(job, config, sample, vcf_files):
    """Use vcf-isec to merge vcfs together and create an ensemble call set report"""

    files = list()
    for vcf in vcf_files:
        utilities.bgzip_and_tabix_vcf(job, vcf)
        files.append("{}.gz".format(vcf))
    vcf_files_string = " ".join(files)

    merged_vcf = "{}.merged.vcf".format(sample)
    logfile = "{}.merged.log".format(sample)

    isec_command = ("{}".format(config['vcftools_isec']['bin']),
                    "-n",
                    "+1",
                    "{}".format(vcf_files_string),
                    ">",
                    "{}".format(merged_vcf))

    job.fileStore.logToMaster("Vcftools intersect Command: {}\n".format(isec_command))
    pipeline.run_and_log_command(" ".join(isec_command), logfile)

    return merged_vcf
