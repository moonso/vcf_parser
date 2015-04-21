import sys

if sys.version_info < (2, 7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict

from vcf_parser.utils import build_info_string

def test_simple_info_string():
    """
    Test how the build_info_string behaves
    """
    
    info = OrderedDict()
    info['MQ'] = ['1']
    
    assert build_info_string(info) == "MQ=1"

def test_info_string_with_no_value():
    """
    Test how the build_info_string behaves
    """
    
    info = OrderedDict()
    info['MQ'] = ['1']
    info['BOOL'] = []
    
    assert build_info_string(info) == "MQ=1;BOOL"
    