from __future__ import unicode_literals
from vcf_parser.utils import build_rank_score_dict
import pytest

def test_build_rank_score_dict():
    """
    Test if the rank score dict is built in proper way.
    """
    
    rank_scores = ['1:12']
    
    assert build_rank_score_dict(rank_scores) == {'1':'12'}


def test_empty_scores():
    """
    Test if the rank score dict is built in proper way.
    """
    
    rank_scores = []
    
    assert build_rank_score_dict(rank_scores) == {}
    
def test_build_multiple_family_rank_score_dict():
    """
    Test if the rank score dict is built in proper way.
    """
    
    rank_scores = ['1:12', '3:4']
    
    assert build_rank_score_dict(rank_scores) == {
        '1':'12', '3':'4'}

def test_malformed_input():
    """
    Test if the rank score dict is built in proper way.
    """
    
    rank_scores = ['1,12', '3:4']
    
    with pytest.raises(SyntaxError):
        build_rank_score_dict(rank_scores)