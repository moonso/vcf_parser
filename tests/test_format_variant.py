from vcf_parser.utils import format_variant
from vcf_parser import Genotype, HeaderParser

import sys
import pytest

if sys.version_info < (2, 7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict

def get_header(header_lines = None):
    """Initiate a HeaderParser and return it"""
    header_parser = HeaderParser()
    
    if not header_lines:
        header_lines = [
            '##fileformat=VCFv4.2',
            '##FILTER=<ID=LowQual,Description="Low quality">',
            '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">',
            '##INFO=<ID=SQ,Number=G,Type=Float,Description="Just for test">',
            '##INFO=<ID=CNT,Number=A,Type=Integer,Description="Number of times '\
            'this allele was found in external db">',
            '##contig=<ID=1,length=249250621,assembly=b37>',
            '##INFO=<ID=DP_HIST,Number=R,Type=String,Description="Histogram for '\
            'DP; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|'\
            '62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5">',
            '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for'\
            ' the ref and alt alleles in the order listed">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=GQ,Number=1,Type=String,Description="GenotypeQuality">'
            '##reference=file:///human_g1k_v37.fasta',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfather\tmother\tproband'
        ]
    for line in header_lines:
        if line.startswith('##'):
            header_parser.parse_meta_data(line)
        elif line.startswith('#'):
            header_parser.parse_header_line(line)
    
    return header_parser


def test_simple_variant():
    """
    Test how the format_variant behaves
    """
    
    header_parser = get_header()
    
    variant_line = "1\t11900\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t"\
                    "0/1:60\t1/1:60"
    
    variant = format_variant(variant_line, header_parser)
    info_dict = OrderedDict()
    info_dict['MQ'] = ['1']
    
    assert variant['CHROM'] == "1"
    assert variant['POS'] == "11900"
    assert variant['ID'] == "."
    assert variant['REF'] == "A"
    assert variant['ALT'] == "T"
    assert variant['QUAL'] == "100"
    assert variant['FILTER'] == "PASS"
    assert variant['INFO'] == "MQ=1"
    assert variant['FORMAT'] == "GT:GQ"
    assert variant['father'] == "0/1:60"
    assert variant['mother'] == "0/1:60"
    assert variant['proband'] == "1/1:60"
    assert variant['info_dict'] == info_dict
    assert type(variant['genotypes']['mother']) == type(Genotype())

def test_malformed_line():
    """
    Test if proper behaviour with malformed vcf line
    """
    
    header_parser = get_header()
    # Missing position
    variant_line = "1\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t"\
                    "0/1:60\t1/1:60"
    
    with pytest.raises(SyntaxError):
        format_variant(variant_line, header_parser)

def test_no_genotypes():
    """
    Test if proper behaviour with minimal vcf
    """
    header_lines = [
        '##fileformat=VCFv4.2',
        '##FILTER=<ID=LowQual,Description="Low quality">',
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">',
        '##contig=<ID=1,length=249250621,assembly=b37>',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
    ]
    
    header_parser = get_header(header_lines)
    # Missing position
    variant_line = "1\t11900\t.\tA\tT\t100\tPASS\tMQ=1"
    
    format_variant(variant_line, header_parser)

    
def test_wrong_number_annotation_integer():
    """
    Test if proper behaviour with malformed vcf line
    """
    
    header_parser = get_header()
    # Missing position
    variant_line = "1\t11900\t.\tA\tT\t100\tPASS\tMQ=1,1\tGT:GQ\t0/1:60\t"\
                    "0/1:60\t1/1:60"
                   
    
    with pytest.raises(SyntaxError):
        format_variant(variant_line, header_parser)
    
def test_wrong_number_annotation_allele():
    """
    Test if proper behaviour with malformed vcf line
    """
    
    header_parser = get_header()
    # Missing position
    variant_line = "3\t947379\t.\tA\tT,C\t100\tPASS\tMQ=1;CNT=5;"\
    "DP_HIST=12,43,22\tGT:GQ:AD:DP\t1/1:60:0,7,0:12\t0/2:60:7,0,10:17"\
    "\t1/2:60:0,7,8:16"
    
    
    with pytest.raises(SyntaxError):
        format_variant(variant_line, header_parser)

def test_wrong_number_annotation_reference():
    """
    Test if proper behaviour with malformed vcf line
    """
    
    header_parser = get_header()
    # Missing position
    variant_line = "3\t947379\t.\tA\tT,C\t100\tPASS\tMQ=1;CNT=5,6;"\
    "DP_HIST=1,2\tGT:GQ:AD:DP\t1/1:60:0,7,0:12\t0/2:60:7,0,10:17"\
    "\t1/2:60:0,7,8:16"
    
    
    with pytest.raises(SyntaxError):
        variant = format_variant(variant_line, header_parser)

def test_wrong_number_annotation_genotype():
    """
    Test if proper behaviour with malformed vcf line
    """
    
    header_parser = get_header()
    # Missing position
    variant_line = "3\t947379\t.\tA\tT,C\t100\tPASS\tMQ=1;CNT=5,6;"\
    "SQ=1,2\tGT:GQ:AD:DP\t1/1:60:0,7,0:12\t0/2:60:7,0,10:17"\
    "\t1/2:60:0,7,8:16"
    
    
    with pytest.raises(SyntaxError):
        variant = format_variant(variant_line, header_parser)
    
    
    
    
    
    
    
    
    
    
    
    