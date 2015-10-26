from vcf_parser.utils import build_vep_string
import pytest

def test_build_vep_string():
    """
    Test how the build_vep_string behaves
    """
    vep_headers = ["Position", "Type", "Alternative"]
    annotation = [
        {"Position":"1", "Type":"nonsynonymous", "Alternative":"C"}
    ]
    
    assert build_vep_string(annotation, vep_headers) == "1|nonsynonymous|C"

def test_null_string():
    """
    Test how the build_vep_string behaves if no annotations
    """
    vep_headers = ["Position", "Type", "Alternative"]
    annotation = []
    
    assert build_vep_string(annotation, vep_headers) == ""

def test_build_multiple_string():
    """
    Test how the build_vep_string behaves
    """
    vep_headers = ["Position", "Type", "Alternative"]
    annotation = [
        {"Position":"1", "Type":"nonsynonymous", "Alternative":"C"},
        {"Position":"2", "Type":"synonymous", "Alternative":"G"},
    ]
    
    assert build_vep_string(annotation, vep_headers) == "1|nonsynonymous|C,"\
                                                        "2|synonymous|G"

def test_wrong_string():
    """
    Test how the build_vep_string behaves
    """
    vep_headers = ["Position", "Type", "Alternative", "Region"]
    annotation = [
        {"Position":"1", "Type":"nonsynonymous", "Alternative":"C"}
    ]
    
    with pytest.raises(SyntaxError):
        build_vep_string(annotation, vep_headers)
