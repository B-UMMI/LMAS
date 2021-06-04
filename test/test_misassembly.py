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
    assert isinstance(paf_dict, dict)
    assert len(paf_dict.keys()) == 12


def test_filter_dict():
    """
    test module to filter dictionary to contain only contigs broken into multiple
    alignment blocks.
    """
    paf_dict = misassembly.parse_paf(MISASSEMBLY_PAF_FILE)
    filter_paf_dict = misassembly.filter_dict(paf_dict)

    assert isinstance(filter_paf_dict, dict)
    assert len(filter_paf_dict.keys()) == 6


def test_classify_misassembled_contigs():
    paf_dict = misassembly.parse_paf(MISASSEMBLY_PAF_FILE)
    filter_paf_dict = misassembly.filter_dict(paf_dict)
    classified_mis_dict = misassembly.classify_misassembled_contigs(
        filter_paf_dict)

    assert sorted(list(classified_mis_dict.keys())) == [
        '162', 'Contig_1255_11.5951', 'NODE_188_length_33202_cov_63.293119','contig-100_82', 'k141_878', 'scaffold_8']

    # classify chimera
    assert 'chimera' in classified_mis_dict['k141_878']['misassembly'][0]
    assert 'Salmonella_enterica' in classified_mis_dict['k141_878']['misassembly'][0]
    assert 'Escherichia_coli_plasmid' in classified_mis_dict['k141_878']['misassembly'][0]
    assert 'Escherichia_coli' in classified_mis_dict['k141_878']['misassembly'][0]

    assert 'chimera' in classified_mis_dict['Contig_1255_11.5951']['misassembly'][0]
    assert 'Pseudomonas_aeruginosa' in classified_mis_dict['Contig_1255_11.5951']['misassembly'][0]
    assert 'Salmonella_enterica' in classified_mis_dict['Contig_1255_11.5951']['misassembly'][0]

    # classify inversion
    assert sorted(classified_mis_dict['contig-100_82']['misassembly']) == [
        'inversion']
    assert len(classified_mis_dict['contig-100_82']['strands']) == 2
    print(classified_mis_dict['contig-100_82'])

    # classify translocation
    assert sorted(classified_mis_dict['scaffold_8']['misassembly']) == [
        'translocation']
    assert len(classified_mis_dict['scaffold_8']['strands']) == 1
    assert any(
        i > 1000 for i in classified_mis_dict['scaffold_8']['distance_in_ref'])

    # classify inversion + translocation
    assert sorted(classified_mis_dict['162']['misassembly']) == [
        'inversion', 'translocation']
    assert len(classified_mis_dict['162']['strands']) == 2
    assert any(i > 1000 for i in classified_mis_dict['162']['distance_in_ref'])

    assert sorted(classified_mis_dict['NODE_188_length_33202_cov_63.293119']['misassembly']) == [
        'inversion', 'translocation']
    assert len(
        classified_mis_dict['NODE_188_length_33202_cov_63.293119']['strands']) == 2
    assert any(
        i > 1000 for i in classified_mis_dict['NODE_188_length_33202_cov_63.293119']['distance_in_ref'])


"""
def test_all():
    paf_dict = misassembly.parse_paf("data/misassembly/EHS_SKESA.paf")
    filter_paf_dict = misassembly.filter_dict(paf_dict)
    classified_mis_dict = misassembly.classify_misassembled_contigs(filter_paf_dict)

    print(len(classified_mis_dict))

    with open("misassembly_test.json", "w") as fh:
        fh.write(str(classified_mis_dict))
"""
