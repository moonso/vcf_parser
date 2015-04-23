from vcf_parser.utils import split_genotype

def test_simple_split():
    """
    Test how split_genotype behaves
    """
    
    genotype = "1/2:30"
    gt_format = "GT:DP"
    
    first_genotype = split_genotype(genotype, gt_format, 0)
    second_genotype = split_genotype(genotype, gt_format, 1)
    third_genotype = split_genotype(genotype, gt_format, 2)
    
    assert first_genotype == "0/1:30"
    assert second_genotype == "0/1:30"
    assert third_genotype == "0/0:30"

def test_simple_split_ad():
    """
    Test how split_genotype behaves with allele depth
    """
    
    genotype = "1/2:30:0,10,20"
    gt_format = "GT:DP:AD"
    
    first_genotype = split_genotype(genotype, gt_format, 0)
    second_genotype = split_genotype(genotype, gt_format, 1)
    third_genotype = split_genotype(genotype, gt_format, 2)
    
    assert first_genotype == "0/1:30:0,10"
    assert second_genotype == "0/1:30:0,20"
    assert third_genotype == "0/0:30:0,0"


def test_other_allele_symbol():
    """
    Test how split_genotype behaves
    """
    
    genotype = "1/2:30"
    gt_format = "GT:DP"
    
    first_genotype = split_genotype(genotype, gt_format, 0, '.')
    second_genotype = split_genotype(genotype, gt_format, 1, '.')
    third_genotype = split_genotype(genotype, gt_format, 2, '.')
    
    assert first_genotype == "./1:30"
    assert second_genotype == "./1:30"
    assert third_genotype == "./.:30"
    