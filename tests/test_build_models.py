from vcf_parser.utils import build_models_dict

def test_simple_models():
    """
    Test how the build_models_dict behaves
    """
    models = ['1:AD','2:AR_comp|AD_dn']
    
    assert build_models_dict(models) == {'1':['AD'], '2':['AR_comp', 'AD_dn']}

def test_empty_models():
    """
    Test how the build_models_dict behaves
    """
    models = []
    
    assert build_models_dict(models) == {}
