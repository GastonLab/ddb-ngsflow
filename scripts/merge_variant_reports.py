#!/usr/bin/env python

import os
import csv
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plot

from collections import defaultdict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samples_file', help="CSV file containing sample meta-data")
    parser.add_argument('-o', '--output', help="Name for summarized output file")
    args = parser.parse_args()

    sample_variants = dict()

    files = os.listdir(os.getcwd())
    for report in files:
        sections = report.partition(".")
        if sections[2] == 'variant_report.txt':
            sample = sections[0]
            with open(report, 'r') as report_file:
                reader = csv.DictReader(report_file, delimiter='\t')
                sample_variants[sample] = list()
                for row in reader:
                    temp_dict = row
                    temp_dict['sample'] = sample
                    sample_variants[sample].append(temp_dict)

    with open(args.output, 'w') as outfile:
        files = os.listdir(os.getcwd())
        fieldnames = ["sample", "chrom", "start", "end", "ref", "alt", "vcf_id", "rs_ids", "cosmic_ids", "filter",
                      "qual", "qual_depth", "depth", "gene", "transcript", "exon", "codon_change", "aa_change",
                      "biotype", "impact", "impact_so", "impact_severity", "aa_length", "is_lof", "is_conserved",
                      "pfam_domain", "in_omim", "clinvar_sig", "clinvar_disease_name", "clinvar_origin",
                      "clinvar_causal_allele", "clinvar_dbsource", "clinvar_dbsource_id", "clinvar_on_diag_assay",
                      "rmsk", "in_segdup", "strand_bias", "rms_map_qual", "in_hom_run", "num_mapq_zero",
                      "num_reads_w_dels", "grc", "gms_illumina", "in_cse", "num_alleles", "allele_count",
                      "haplotype_score", "is_somatic", "somatic_score", "aaf_esp_ea", "aaf_esp_aa", "aaf_esp_aa",
                      "aaf_esp_all", "aaf_1kg_amr", "aaf_1kg_eas", "aaf_1kg_sas", "aaf_1kg_afr", "aaf_1kg_eur",
                      "aaf_1kg_all", "aaf_exac_all", "aaf_adj_exac_all", "aaf_adj_exac_afr", "aaf_adj_exac_amr",
                      "aaf_adj_exac_eas", "aaf_adj_exac_fin", "aaf_adj_exac_nfe", "aaf_adj_exac_oth",
                      "aaf_adj_exac_sas", "max_aaf_all", "in_esp", "in_1kg", "in_exac"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for sample in sample_variants:
            for row in sample_variants[sample]:
                writer.writerow(row)

