#!/usr/bin/env python

import sys
import cyvcf
from cyvcf import utils

vcf_readers = list()
for vcffile in sys.argv[1:]:
    vcf_reader = cyvcf.Reader(open(vcffile, 'rb'))
    vcf_readers.append(vcf_reader)

for caller_records in utils.walk_together(*vcf_readers):
    num_callers = 0
    for record in caller_records:
        if record is not None:
            num_callers += 1
    if num_callers > 1:
        for record in caller_records:
            sys.stdout.write("*")
            print record
        sys.stdout.write("\n")
