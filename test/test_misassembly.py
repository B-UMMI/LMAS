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
MISASSEMBLY_PAF_FILE_ALL = "test/data/misassembly/test.paf"
MISASSEMBLY_PAF_FILE_CHIMERA = "test/data/misassembly/test_chimera.paf"
MISASSEMBLY_PAF_FILE_INVERSION = "test/data/misassembly/test_inversion.paf"
MISASSEMBLY_PAF_FILE_INSERTION = "test/data/misassembly/test_insertion.paf"
MISASSEMBLY_PAF_FILE_TRANSLOCATION = "test/data/misassembly/test_translocation.paf"
MISASSEMBLY_PAF_FILE_COMPLEX = "test/data/misassembly/test_inversion_translocation.paf"


def test_parse_paf():
    """
    test module to parse paf files
    """
    paf_dict = misassembly.parse_paf(MISASSEMBLY_PAF_FILE_ALL)
    assert isinstance(paf_dict, dict)
    assert len(paf_dict.keys()) == 13


def test_filter_dict():
    """
    test module to filter dictionary to contain only contigs broken into multiple
    alignment blocks.
    """
    paf_dict = misassembly.parse_paf(MISASSEMBLY_PAF_FILE_ALL)
    filter_paf_dict = misassembly.filter_dict(paf_dict)

    assert isinstance(filter_paf_dict, dict)
    assert len(filter_paf_dict.keys()) == 7


def test_chimera():
    """
    Test the detection of chimeric contigs
    """
    paf_dict = misassembly.parse_paf(MISASSEMBLY_PAF_FILE_CHIMERA)
    filter_paf_dict = misassembly.filter_dict(paf_dict)
    classified_mis_dict = misassembly.classify_misassembled_contigs(
        filter_paf_dict)

    assert sorted(list(classified_mis_dict.keys())) == [
        'Contig_1255_11.5951', 'k141_878', ]

    # classify chimera
    assert 'chimera' in classified_mis_dict['k141_878']['misassembly'][0]
    assert 'Salmonella_enterica' in classified_mis_dict['k141_878']['misassembly'][0]
    assert 'Escherichia_coli_plasmid' in classified_mis_dict['k141_878']['misassembly'][0]
    assert 'Escherichia_coli' in classified_mis_dict['k141_878']['misassembly'][0]

    assert 'chimera' in classified_mis_dict['Contig_1255_11.5951']['misassembly'][0]
    assert 'Pseudomonas_aeruginosa' in classified_mis_dict['Contig_1255_11.5951']['misassembly'][0]
    assert 'Salmonella_enterica' in classified_mis_dict['Contig_1255_11.5951']['misassembly'][0]


def test_inversion():
    """
    Test the detection of inversions in the contigs
    """
    paf_dict = misassembly.parse_paf(MISASSEMBLY_PAF_FILE_INVERSION)
    filter_paf_dict = misassembly.filter_dict(paf_dict)
    classified_mis_dict = misassembly.classify_misassembled_contigs(
        filter_paf_dict)

    assert sorted(list(classified_mis_dict.keys())) == [
        'NODE_84_length_7746_cov_1230.235860', 'NODE_85_length_7744_cov_815.448823', 'contig-100_82']

    for contig in classified_mis_dict.keys():
        assert len(list(classified_mis_dict[contig]['strands'])) > 1


def test_insertion():
    """
    Test the detection of insertions in the contigs (sequences in contigs that are not present/don't map to the reference)
    """
    paf_dict = misassembly.parse_paf(MISASSEMBLY_PAF_FILE_INSERTION)
    filter_paf_dict = misassembly.filter_dict(paf_dict)
    classified_mis_dict = misassembly.classify_misassembled_contigs(
        filter_paf_dict)
    
    assert sorted(list(classified_mis_dict.keys())) == [
        'NODE_55_length_174716_cov_35.030820', 'scaffold_6']

    assert sorted(classified_mis_dict['NODE_55_length_174716_cov_35.030820']['misassembly']) == [
        'insertion', 'rearrangement']
    assert any(
        i > 50 for i in classified_mis_dict['NODE_55_length_174716_cov_35.030820']['distance_in_contig'])

    assert sorted(classified_mis_dict['scaffold_6']['misassembly']) == [
        'insertion']
    assert any(
        i > 50 for i in classified_mis_dict['scaffold_6']['distance_in_contig'])


def test_translocation():
    "test the detection of translocation"

    paf_dict = misassembly.parse_paf(MISASSEMBLY_PAF_FILE_TRANSLOCATION)
    filter_paf_dict = misassembly.filter_dict(paf_dict)
    classified_mis_dict = misassembly.classify_misassembled_contigs(
        filter_paf_dict)
    
    # classify translocation
    assert sorted(classified_mis_dict['scaffold_8']['misassembly']) == [
        'insertion', 'rearrangement', 'translocation']
    assert len(classified_mis_dict['scaffold_8']['strands']) == 1
    assert any(
        i > 1000 for i in classified_mis_dict['scaffold_8']['distance_in_ref'])


def test_complex():
    paf_dict = misassembly.parse_paf(MISASSEMBLY_PAF_FILE_COMPLEX)
    filter_paf_dict = misassembly.filter_dict(paf_dict)
    classified_mis_dict = misassembly.classify_misassembled_contigs(
        filter_paf_dict)
    
    assert sorted(classified_mis_dict['162']['misassembly']) == [
        'inversion', 'rearrangement', 'translocation']
    assert len(classified_mis_dict['162']['strands']) == 2
    assert any(i > 1000 for i in classified_mis_dict['162']['distance_in_ref'])

    assert sorted(classified_mis_dict['NODE_188_length_33202_cov_63.293119']['misassembly']) == [
        'inversion', 'rearrangement', 'translocation']
    assert len(
        classified_mis_dict['NODE_188_length_33202_cov_63.293119']['strands']) == 2
    assert any(
        i > 1000 for i in classified_mis_dict['NODE_188_length_33202_cov_63.293119']['distance_in_ref'])

"""
def test_make_df():
    paf_dict = misassembly.parse_paf(MISASSEMBLY_PAF_FILE_ALL)
    filter_paf_dict = misassembly.filter_dict(paf_dict)
    classified_mis_dict = misassembly.classify_misassembled_contigs(
        filter_paf_dict)

    print(classified_mis_dict)
    reference_report = {"sample": "lala",
                        "assembler": "lala", 'reference': {}}
    for contig in classified_mis_dict.keys():
        print(classified_mis_dict[contig].keys())
        for reference in classified_mis_dict[contig]['reference']:
            if reference not in reference_report['reference'].keys():
                reference_report['reference'][reference] = 1
            else:
                reference_report['reference'][reference] += 1

    
    #print(filter_paf_dict)
    #print(classified_mis_dict)

    #for contig_id in classified_mis_dict.keys():
    #    for contig_info in filter_paf_dict[contig_id]:
    #        print(contig_info)
    
    #bubu = misassembly.make_df('pytest_sample', 'pytest', classified_mis_dict, filter_paf_dict)
    #print(bubu)
"""
