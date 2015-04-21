from vcf_parser.utils import format_variant
from vcf_parser import Genotype

import sys

if sys.version_info < (2, 7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict

def test_simple_variant():
    """
    Test how the format_variant behaves
    """
    variant_line = "1\t11900\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t"\
                    "0/1:60\t1/1:60"
    individuals = ['father','mother','proband']
    vcf_header = [
        'CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',
        'father','mother','proband'
    ]
    vep_header = []
    variant = format_variant(variant_line, individuals, vcf_header, vep_header)
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

