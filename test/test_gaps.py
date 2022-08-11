#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Module to test the LMAS script to detect and classify misassemblies

Raises
------
pytest.fail
    Status for the test (Pass or Fail)
"""
from templates.gap_assessment import get_gaps
import pytest
from contextlib import contextmanager
try:
    from templates import gap_assessment
except ImportError:
    from templates import gap_assessment


@contextmanager
def not_raises(exception, msg):
    """[summary]

    Parameters
    ----------
    exception : [type]
        [description]
    msg : [type]
        [description]

    Raises
    ------
    pytest.fail
        [description]
    """
    try:
        yield
    except exception:
        raise pytest.fail(msg)

GAPS_PAF = "test/data/gaps/tiny.paf"
GAPS_PAF_BAD = "test/data/gaps/tiny_bad.paf"

def test_get_gaps_simple():

    gaps, gap_sizes = gap_assessment.get_gaps(GAPS_PAF,"Bacillus_subtilis", 12137031/3)
    
    for i in range(0,len(gaps)-1):
        assert gaps[i][1] > gaps[i][0]
        assert gaps[i][1] < gaps[i+1][0]
        assert len(gaps) == 173

    assert gaps[0][0] > 0
    assert gaps[len(gaps)-1][1] < 12137031/3

    assert len(gap_sizes) == len(gaps)
    for gap_size in gap_sizes:
        assert gap_size > 0 


def test_get_gaps_bad_assembly():

    gaps, gap_sizes = gap_assessment.get_gaps(GAPS_PAF_BAD,"Bacillus_subtilis", 12137031/3)
    
    for i in range(0,len(gaps)-1):
        assert gaps[i][1] > gaps[i][0]
        assert gaps[i][1] < gaps[i+1][0]
        assert len(gaps) == 60

    assert gaps[0][0] > 0
    assert gaps[len(gaps)-1][1] < 12137031/3

    assert len(gap_sizes) == len(gaps)
    for gap_size in gap_sizes:
        assert gap_size > 0 
