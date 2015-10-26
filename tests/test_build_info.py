import sys

if sys.version_info < (2, 7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict

from vcf_parser.utils import build_info_string, build_info_dict

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


def test_build_simple_info_dict():
    """
    Test to build a simple info dict from a info string
    """
    
    info_string = "MQ=1;CNT=5,8;DP_HIST=12,43,22;RS"
    info_dict = OrderedDict()
    info_dict['MQ'] = ['1']
    info_dict['CNT']= ['5','8']
    info_dict['DP_HIST'] = ['12','43','22']
    info_dict['RS'] = []
    
    assert build_info_dict(info_string) == info_dict