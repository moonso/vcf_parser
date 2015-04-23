#!/usr/bin/env python
# encoding: utf-8
"""
test_variant_parser.py

Test the so that the vcf parser behave as suspected.

Created by Måns Magnusson on 2014-03-04.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.
"""

from tempfile import NamedTemporaryFile
from vcf_parser import VCFParser
import pytest



def get_vcf_file(vcf_lines):
    """
    Take an iterator with vcf lines and prints them to a temporary file.
    
    Arguments:
        vcf_lines (iterator): An iterator with vcf lines
    
    Returns:
        filename (str): The path to the vcf file
    """
    vcf_file = NamedTemporaryFile(mode='w+t', delete=False, suffix='.vcf')
    vcf_file.writelines(vcf_lines)
    vcf_file.seek(0)
    vcf_file.close()
    
    return vcf_file.name
    

def test_vcf_parser():
    """
    Test the vcf_parser
    """
    vcf_lines = [
        '##fileformat=VCFv4.1\n',
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">\n',
        '##contig=<ID=1,length=249250621,assembly=b37>\n',
        '##reference=file:///humgen/gsa-hpprojects/GATK/bundle'\
        '/current/b37/human_g1k_v37.fasta\n',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'\
        'father\tmother\tproband\n',
        '1\t11900\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/1:60\t1/1:60\n',
        '1\t879585\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/0:60\t0/1:60\n',
        '1\t879586\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/1:60\t0/1:60\n',
        '1\t947378\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/0:60\t0/1:60\n',
        '1\t973348\t.\tG\tA\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/0:60\t0/1:60\n',
        '3\t879585\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/0:60\t0/1:60\n',
        '3\t879586\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/1:60\t0/1:60\n',
        '3\t947378\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/0:60\t0/1:60\n',
        '3\t973348\t.\tG\tA\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/0:60\t0/1:60\n'
        ]
    
    vcf_file = get_vcf_file(vcf_lines)
    variants = []
    for variant in VCFParser(vcf_file):
        variants.append(variant)
    
    first_variant = variants[0]
    assert first_variant['POS'] == '11900'
    
    second_variant = variants[1]
    assert second_variant['POS'] == '879585'
    
    last_variant = variants[-1]
    assert last_variant['POS'] == '973348'

def test_one_variant():
    """
    Test the vcf_parser
    """
    vcf_lines = [
        '##fileformat=VCFv4.1\n',
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">\n',
        '##contig=<ID=1,length=249250621,assembly=b37>\n',
        '##reference=file:///humgen/gsa-hpprojects/GATK/bundle'\
        '/current/b37/human_g1k_v37.fasta\n',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'\
        'father\tmother\tproband\n',
        '1\t11900\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/1:60\t1/1:60\n',
        ]
    
    vcf_file = get_vcf_file(vcf_lines)
    variants = []
    for variant in VCFParser(infile=vcf_file):
        variants.append(variant)
    
    first_variant = variants[0]
    assert first_variant['POS'] == '11900'
    assert first_variant['ALT'] == 'T'
    

def test_split_variant():
    """
    Test the vcf_parser
    """
    vcf_lines = [
        '##fileformat=VCFv4.1\n',
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">\n',
        '##contig=<ID=1,length=249250621,assembly=b37>\n',
        '##reference=file:///humgen/gsa-hpprojects/GATK/bundle'\
        '/current/b37/human_g1k_v37.fasta\n',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'\
        'father\tmother\tproband\n',
        '1\t11900\t.\tA\tT,C\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/2:60\t1/2:60\n',
        ]
    
    vcf_file = get_vcf_file(vcf_lines)
    variants = []
    for variant in VCFParser(infile=vcf_file, split_variants=True):
        variants.append(variant)
    
    first_variant = variants[0]
    assert first_variant['POS'] == '11900'
    assert first_variant['ALT'] == 'T'
    
    second_variant = variants[1]
    assert second_variant['POS'] == '11900'
    assert second_variant['ALT'] == 'C'


def test_wrong_formatted_vcf():
    """
    Test how vcf_parser behaves if no fileformat is given
    """
    vcf_lines = [
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">\n',
        '##contig=<ID=1,length=249250621,assembly=b37>\n',
        '##reference=file:///humgen/gsa-hpprojects/GATK/bundle'\
        '/current/b37/human_g1k_v37.fasta\n',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'\
        'father\tmother\tproband\n',
        '1\t11900\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/1:60\t1/1:60\n',
        '1\t879585\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/0:60\t0/1:60\n',
        '1\t879586\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/1:60\t0/1:60\n',
        '1\t947378\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/0:60\t0/1:60\n',
        '1\t973348\t.\tG\tA\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/0:60\t0/1:60\n',
        '3\t879585\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/0:60\t0/1:60\n',
        '3\t879586\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/1:60\t0/1:60\n',
        '3\t947378\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/0:60\t0/1:60\n',
        '3\t973348\t.\tG\tA\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/0:60\t0/1:60\n'
        ]
    vcf_file = get_vcf_file(vcf_lines)
    with pytest.raises(SyntaxError):
        for variant in VCFParser(vcf_file):
            print(variant)
    
def test_build_vcf():
    """
    Test how it works to build a vcf by adding metadata and variants to the parser
    """
    parser = VCFParser()
    variants = []
    parser.metadata.add_fileformat("VCFv4.1")
    assert parser.metadata.fileformat == "VCFv4.1"
    
    parser.metadata.add_info(
        info_id='MQ', number='1', entry_type='Float', description="RMS Mapping Quality")
    
    assert 'MQ' in parser.metadata.extra_info
    
    # No header line added yet
    with pytest.raises(SyntaxError):
        parser.add_variant(chrom='1', pos='11900', rs_id='.', ref='A', 
                        alt='T', qual='100', filt='PASS', info="MQ=1")
    
    header_line = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'\
        'father\tmother\tproband'
    parser.metadata.parse_header_line(header_line)
    parser.add_variant(chrom='1', pos='11900', rs_id='.', ref='A', 
                        alt='T', qual='100', filt='PASS', info="MQ=1",
                        form="GT:GQ", genotypes=["0/1:60", "0/1:60", "1/1:60"])
    for variant in parser:
        variants.append(variant)
    
    first_variant = variants[0]
    
    assert first_variant['POS'] == '11900'

