#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Module to test the LMAS script to detect and classify misassemblies

Raises
------
pytest.fail
    Status for the test (Pass or Fail)
"""
from templates.misassembly import classify_misassembled_contigs
import pytest
from contextlib import contextmanager
try:
    from templates import misassembly
except ImportError:
    from templates import misassembly


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

# GLOBAL VARIABLES - TEST FILES
MISASSEMBLY_PAF_FILE = "data/misassembly/test.paf"


def test_parse_paf():
    """
    test module to parse paf files
    """
    paf_dict = misassembly.parse_paf(MISASSEMBLY_PAF_FILE)
    assert isinstance(paf_dict,dict)
    assert len(paf_dict.keys()) == 9

def test_filter_dict():
    """
    test module to filter dictionary to contain only contigs broken into multiple
    alignment blocks.
    """
    paf_dict = misassembly.parse_paf(MISASSEMBLY_PAF_FILE)
    filter_paf_dict = misassembly.filter_dict(paf_dict)

    assert isinstance(filter_paf_dict,dict)
    assert len(filter_paf_dict.keys()) == 3

def test_classify_misassembled_contigs():
    paf_dict = misassembly.parse_paf(MISASSEMBLY_PAF_FILE)
    filter_paf_dict = misassembly.filter_dict(paf_dict)
    classified_mis_dict = misassembly.classify_misassembled_contigs(filter_paf_dict)

    assert isinstance(classified_mis_dict,dict)
    assert len(classified_mis_dict.keys()) == 3

    assert list(classified_mis_dict.keys()) == ['162', 'NODE_188_length_33202_cov_63.293119',  'k141_878']

    assert classified_mis_dict['162']['misassembly'] == ['inversion', 'translocation']
    assert classified_mis_dict['NODE_188_length_33202_cov_63.293119']['misassembly'] == ['inversion', 'translocation']
    assert 'chimera' in classified_mis_dict['k141_878']['misassembly'][0]
    assert 'Salmonella_enterica' in classified_mis_dict['k141_878']['misassembly'][0]
    assert 'Escherichia_coli_plasmid' in classified_mis_dict['k141_878']['misassembly'][0]
    assert 'Escherichia_coli' in classified_mis_dict['k141_878']['misassembly'][0]


def test_all():
    paf_dict = misassembly.parse_paf("data/misassembly/EHS_GATBMiniaPipeline.paf")
    filter_paf_dict = misassembly.filter_dict(paf_dict)
    classified_mis_dict = misassembly.classify_misassembled_contigs(filter_paf_dict)

    print(len(classified_mis_dict))

    with open("misassembly_test.json", "w") as fh:
        fh.write(str(classified_mis_dict))