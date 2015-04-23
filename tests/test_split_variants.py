from vcf_parser.utils import split_variants, format_variant
from vcf_parser import HeaderParser

def get_header(header_lines = None):
    """Initiate a HeaderParser and return it"""
    header_parser = HeaderParser()
    
    if not header_lines:
        header_lines = [
            '##fileformat=VCFv4.2',
            '##FILTER=<ID=LowQual,Description="Low quality">',
            '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">',
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


def test_simple_split():
    """
    Test how split genotypes behave when a simple split
    """
    
    header_parser = get_header()
    
    variant_line = "3\t947379\t.\tA\tT,C\t100\tPASS\tMQ=1;CNT=5,8;"\
    "DP_HIST=12,43,22\tGT:GQ:AD:DP\t1/1:60:0,7,0:12\t0/2:60:7,0,10:17"\
    "\t1/2:60:0,7,8:16"
    
    variant = format_variant(variant_line, header_parser)
    
    splitted_variants = []
    
    for variant in split_variants(variant, header_parser):
        splitted_variants.append(variant)
    
    assert len(splitted_variants) == 2
     
    first_variant = splitted_variants[0]
    second_variant = splitted_variants[1]
    
    # Test if the splitted variants still have the same reference
    assert first_variant['REF'] == 'A'
    assert second_variant['REF'] == 'A'
    # Test if the alternative was splitted properly
    assert first_variant['ALT'] == 'T'
    assert second_variant['ALT'] == 'C'
    # Test if simple ino field is handled correct
    assert first_variant['info_dict']['MQ'] == ['1']
    assert second_variant['info_dict']['MQ'] == ['1']
    # Test if info field with Number='A' is handled correct
    assert first_variant['info_dict']['CNT'] == ['5']
    assert second_variant['info_dict']['CNT'] == ['8']
    # Test if info field with Number='R' is handled correct
    assert first_variant['info_dict']['DP_HIST'] == ['12', '43']
    assert second_variant['info_dict']['DP_HIST'] == ['12', '22']
    
    # Test if the genortypes are on the correct format
    assert first_variant['father'] == "1/1:60:0,7:12"
    assert second_variant['father'] == "0/0:60:0,0:12"
    
    assert first_variant['mother'] == "0/0:60:7,0:17"
    assert second_variant['mother'] == "0/1:60:7,10:17"
    
    assert first_variant['proband'] == "0/1:60:0,7:16"
    assert second_variant['proband'] == "0/1:60:0,8:16"

def test_split_minimal():
    """
    Test to split a vcf line without genotypes
    """
    header_lines = [
        '##fileformat=VCFv4.2',
        '##FILTER=<ID=LowQual,Description="Low quality">',
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">',
        '##contig=<ID=1,length=249250621,assembly=b37>',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
    ]
    
    header_parser = get_header(header_lines)
    
    variant_line = "3\t947379\t.\tA\tT,C\t100\tPASS\tMQ=1"
    
    variant = format_variant(variant_line, header_parser)
    
    splitted_variants = []
    
    for variant in split_variants(variant, header_parser):
        splitted_variants.append(variant)
    
    assert len(splitted_variants) == 2
    