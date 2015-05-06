from vcf_parser.utils import build_vep_annotation
import pytest

def test_build_vep_dict():
    """
    Test how the build_vep_dict behaves
    """
    vep_headers = ["Position", "Type", "Allele"]
    annotation = ["1|nonsynonymous|C"]
    
    vep_dict = build_vep_annotation(
        csq_info=annotation, 
        reference='A', 
        alternatives=['C'], 
        vep_columns=vep_headers
        )

    assert vep_dict['C'] == [{
                                "Position":"1",
                                "Type":"nonsynonymous",
                                "Allele":"C"
                            }]

    assert vep_dict['gene_ids'] == set(['-'])


def test_wrong_formated_vep_annotation():
    """
    Test how the build_vep_dict behaves
    """
    vep_headers = ["Position", "Type", "Allele"]
    annotation = ["1|nonsynonymous"]
    
    with pytest.raises(SyntaxError):
        vep_dict = build_vep_annotation(
            csq_info=annotation, 
            reference='A', 
            alternatives=['C'], 
            vep_columns=vep_headers
            )

def test_build_multiple_string():
    """
    Test how the build_vep_string behaves
    """
    vep_headers = ["Position", "Type", "Allele"]
    annotation = [
        "1|nonsynonymous|C",
        "1|synonymous|G",
    ]
    
    vep_dict = build_vep_annotation(
        csq_info=annotation, 
        reference='A', 
        alternatives=['C','G'], 
        vep_columns=vep_headers
        )
    
    assert vep_dict['C'] == [{
                                "Position":"1",
                                "Type":"nonsynonymous",
                                "Allele":"C"
                            }]

    assert vep_dict['G'] == [{
                                "Position":"1",
                                "Type":"synonymous",
                                "Allele":"G"
                            }]

def test_gene_ids():
    """
    Test how the build_vep_string behaves
    """
    vep_headers = ["Position", "Type", "Allele", "SYMBOL"]
    annotation = [
        "1|nonsynonymous|C|ADK",
        "1|synonymous|G|ADK",
    ]
    
    vep_dict = build_vep_annotation(
        csq_info=annotation, 
        reference='A', 
        alternatives=['C','G'], 
        vep_columns=vep_headers
        )
    
    assert vep_dict['gene_ids'] == set(["ADK"])

def test_insertion():
    """
    Test how build vep string behaves with a insertion
    """
    vep_headers = ["Allele", "Feature_type", "Consequence", "SYMBOL"]
    annotation = [
        "ACC|Transcript|downstream_gene_variant|NOC2L", 
        ]
    
    vep_dict = build_vep_annotation(
        csq_info=annotation, 
        reference='C', 
        alternatives=['CACC'], 
        vep_columns=vep_headers
        )
    
    assert vep_dict['CACC'] == [
        {'Allele': 'ACC',
         'Consequence': 'downstream_gene_variant',
         'Feature_type': 'Transcript',
         'SYMBOL': 'NOC2L'
         }
      ]


def test_deletion():
    """
    Test how build vep string behaves with a deletion
    """
    vep_headers = ["Allele", "Feature_type", "Consequence", "SYMBOL"]
    annotation = [
        "-|Transcript|downstream_gene_variant|NOC2L", 
        "-|Transcript|inframe_deletion|SAMD11"
        ]
    
    vep_dict = build_vep_annotation(
        csq_info=annotation, 
        reference='ACC', 
        alternatives=['C'], 
        vep_columns=vep_headers
        )
    
    assert vep_dict['C'] == [
        {'Allele': '-',
         'Consequence': 'downstream_gene_variant',
         'Feature_type': 'Transcript',
         'SYMBOL': 'NOC2L'
         },
         {'Allele': '-',
          'Consequence': 'inframe_deletion',
          'Feature_type': 'Transcript',
          'SYMBOL': 'SAMD11'
          }
      ]


def test_deletion_not_all_bases():
    """
    Test how build vep string behaves with a deletion
    """
    vep_headers = ["Allele", "Feature_type", "Consequence", "SYMBOL"]
    annotation = [
        "TCA|Transcript|inframe_deletion|NOC2L", 
        "TCA|Transcript|downstream_gene_variant|SAMD11"
        ]
    
    vep_dict = build_vep_annotation(
        csq_info=annotation, 
        reference='CTCATCA', 
        alternatives=['CTCA'], 
        vep_columns=vep_headers
        )
    
    assert vep_dict['CTCA'] == [
        {'Allele': 'TCA',
         'Consequence': 'inframe_deletion',
         'Feature_type': 'Transcript',
         'SYMBOL': 'NOC2L'
         },
         {'Allele': 'TCA',
          'Consequence': 'downstream_gene_variant',
          'Feature_type': 'Transcript',
          'SYMBOL': 'SAMD11'
          }
      ]

def test_one_deletion_and_substitution():
    """
    Test how build vep string behaves with a deletion and substitution on
    same position
    """
    vep_headers = ["Allele", "Feature_type", "SYMBOL"]
    annotation = [
        "ACT|Transcript|NOC2L", 
        "-|Transcript|NOC2L"
        ]
    
    vep_dict = build_vep_annotation(
        csq_info=annotation, 
        reference='ACC', 
        alternatives=['ACT', 'C'], 
        vep_columns=vep_headers
        )
    
    assert vep_dict['ACT'] == [
        {'Allele': 'ACT',
         'Feature_type': 'Transcript',
         'SYMBOL': 'NOC2L'
         }]
         
    assert vep_dict['C'] == [
        {'Allele': '-',
         'Feature_type': 'Transcript',
         'SYMBOL': 'NOC2L'
          }
      ]

def test_multiple_deletions():
    """
    Test how build vep string behaves with multiple deletions
    """
    vep_headers = ["Allele", "Feature_type", "SYMBOL"]
    annotation = [
        "ACTT|Transcript|NOC2L", 
        "ACT|Transcript|NOC2L"
        ]
    
    vep_dict = build_vep_annotation(
        csq_info=annotation, 
        reference='TACTTT', 
        alternatives=['TACT', 'TACTT'], 
        vep_columns=vep_headers
        )
    
    assert vep_dict['TACT'] == [
        {'Allele': 'ACT',
         'Feature_type': 'Transcript',
         'SYMBOL': 'NOC2L'
         }]
         
    assert vep_dict['TACTT'] == [
        {'Allele': 'ACTT',
         'Feature_type': 'Transcript',
         'SYMBOL': 'NOC2L'
          }
      ]
