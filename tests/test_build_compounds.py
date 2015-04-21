from vcf_parser.utils import build_compounds_dict

def test_build_compounds():
    """
    Test how build_compounds_dict behaves
    """
    
    compound_info = ['2:1_24_T_A']
    
    assert build_compounds_dict(compound_info) == {
        '2':[
            {
                'variant_id': '1_24_T_A',
                'compound_score': None
            }
        ]
    }

def test_build_multiple_compounds():
    """
    Test how build_compounds_dict behaves
    """
    
    compound_info = ['2:1_24_T_A|1_25_A_C']
    
    assert build_compounds_dict(compound_info) == {
        '2':[
            {
                'variant_id': '1_24_T_A',
                'compound_score': None
            },
            {
                'variant_id': '1_25_A_C',
                'compound_score': None
            }
        ]
    }

def test_build_compounds_with_score():
    """
    Test how build_compounds_dict behaves
    """
    
    compound_info = ['2:1_24_T_A>17']
    
    assert build_compounds_dict(compound_info) == {
        '2':[
            {
                'variant_id': '1_24_T_A',
                'compound_score': '17'
            }
        ]
    }
    