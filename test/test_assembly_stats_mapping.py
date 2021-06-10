#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Module to test the LMAS script to get mapping stats from an assembly file (FASTA)
and a mapping file (PAF)

Raises
------
pytest.fail
    Status for the test (Pass or Fail)
"""
import csv
import pytest
from contextlib import contextmanager
from templates import assembly_stats_global
from templates import utils

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
ASSEMBLY_TEST = "test/data/assembly.fasta"
MAPPING_TEST = "test/data/assembly.paf"

def test_get_mapped_contigs_with_ref():

    mapped_contigs = utils.get_mapped_contigs_with_ref(MAPPING_TEST)
    assert len(mapped_contigs) == 6704

def test_parse_assemblies():

    df = utils.parse_assemblies("pytest", "pytest", ASSEMBLY_TEST, MAPPING_TEST)
    assert df.shape == (6704, 7)
    print(df.columns)