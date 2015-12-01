__author__ = 'dgaston'

import argparse
import pyexcel
import pyexcel.ext.xlsx

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', help="Coverage report")
    parser.add_argument('-o', '--outfile', help="Output file name")
    args = parser.parse_args()

    sheet = pyexcel.get_records(file_name=args.infile, name_columns_by_row=0)
    records = sheet.to_records()
    for record in records:
        low_coverage_samples = list()
        coverage_gaps_samples = list()
        keys = sorted(record.keys())
        for key in keys:
            if 'diagnosetargets.vcf_filter' in key:
                if 'LOW_COVERAGE' in record[key]:
                    low_coverage_samples.append(key)
                if 'COVERAGE_GAPS' in record[key]:
                    coverage_gaps_samples.append(key)
