#!/usr/bin/env python
# encoding: utf-8
"""
test_variant_parser.py

Test the so that the vcf parser behave as suspected.

Created by Måns Magnusson on 2014-03-04.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from tempfile import NamedTemporaryFile
from pprint import pprint as pp

from vcf_parser import VCFParser

class TestVariantParser(object):
    """Test how how the variant parser behaves."""

    def setup_class(self):
        """Setup a vcf file with some variants and a interval tree with features."""
        vcf_lines = [
            '##fileformat=VCFv4.1\n', 
            '##contig=<ID=1,length=249250621,assembly=b37>\n',
            '##reference=file:///humgen/gsa-hpprojects/GATK/bundle'\
            '/current/b37/human_g1k_v37.fasta\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'\
            'father\tmother\tproband\n',
            '1\t11900\t.\tA\tT\t100\tPASS\t MQ=1\t GT:GQ\t0/1:60\t0/1:60\t1/1:60\n',
            '1\t879585\t.\tA\tT\t100\tPASS\t MQ=1\t GT:GQ\t0/1:60\t0/0:60\t0/1:60\n',
            '1\t879586\t.\tA\tT\t100\tPASS\t MQ=1\t GT:GQ\t0/0:60\t0/1:60\t0/1:60\n',
            '1\t947378\t.\tA\tT\t100\tPASS\t MQ=1\t GT:GQ\t0/0:60\t0/0:60\t0/1:60\n',
            '1\t973348\t.\tG\tA\t100\tPASS\t MQ=1\t GT:GQ\t0/0:60\t0/0:60\t0/1:60\n',
            '3\t879585\t.\tA\tT\t100\tPASS\t MQ=1\t GT:GQ\t0/1:60\t0/0:60\t0/1:60\n',
            '3\t879586\t.\tA\tT\t100\tPASS\t MQ=1\t GT:GQ\t0/0:60\t0/1:60\t0/1:60\n',
            '3\t947378\t.\tA\tT\t100\tPASS\t MQ=1\t GT:GQ\t0/0:60\t0/0:60\t0/1:60\n',
            '3\t973348\t.\tG\tA\t100\tPASS\t MQ=1\t GT:GQ\t0/0:60\t0/0:60\t0/1:60\n'
            ]
        self.vcf_file = NamedTemporaryFile(mode='w+t', delete=False, suffix='.vcf')
        self.vcf_file.writelines(vcf_lines)
        self.vcf_file.seek(0)
        self.vcf_file.close()
        
    def test_parser(self):
        """docstring for test_hej"""
        for variant in VCFParser(self.vcf_file.name):
            print(variant)

def main():
    pass


if __name__ == '__main__':
    main()

