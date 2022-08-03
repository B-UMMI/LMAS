#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Module to test the LMAS script to detect SNPs

Raises
------
pytest.fail
    Status for the test (Pass or Fail)
"""
import pytest
import re
from contextlib import contextmanager

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
CIGAR_PERFECT_MATCH = "cs:Z::8029"
CIGAR_SUBSTITUTION = "cs:Z::935*ct:32*ct:117*ga:23*at:204"
CIGAR_SUBSTITUTION_DELETION = "cs:Z::329-gttgctccggttgctcccgttactccagttactccagttgcgccagcttccccggttgcaccagttactcctgttgctccagttgatcctgttaaccctgttgatcctgtttctcccgttgctccagttggacctgttgatcccgttacgccagttgctccggttgcaccagttgttcccgttgctcctgttgatcccgttgctccagtttccccggttgcaccagttgatcccgttgctcctgttgatcctgtttccccggttgcaccagttactcccgttgctccagttgcaccagtggcaccagttactcccgttactcccgttactcccgttactcccgttgcgccagtttccccggttgcaccagttgatccc:5*ta:2*ta:4*at:17*ca:3*ac:25184"


def get_snps(cigar):
    """
    Iterator to get snps from a paf file
    Obtained from templates/snp_assessment.py get_snps() and get_position() methods
    """
    if len(re.findall(r'\*', cigar)) > 0:
        matches = re.findall(r'([:=*+-])(\\d+|[A-Za-z]+)', cigar)
        coord = 0
        substitution = []
        for match in matches:
            if match[0] == ':':
                if match[1] == 'Z':
                    continue
                if int(match[1]): #exact matches to the ref
                    coord += int(match[1])
            if match[0] == '-': # if deletion, coords in the ref need to be ajusted
                coord += len(match[1])
            #if match[0] == '+': # insertions to the reference don't matter
                #coord += len(match[1])
            if match[0] == '*': 
                substitution = match[1]
                coord += 1
                yield coord, substitution

def test_perfect_match():
    snps = get_snps(CIGAR_PERFECT_MATCH)
    assert sum(1 for _ in snps) == len(re.findall(r'\*', CIGAR_PERFECT_MATCH))

def test_substitution():
    snps = get_snps(CIGAR_SUBSTITUTION)
    assert sum(1 for _ in snps) == len(re.findall(r'\*', CIGAR_SUBSTITUTION))

def test_substitution_deletion():
    snps = get_snps(CIGAR_SUBSTITUTION_DELETION)
    assert sum(1 for _ in snps) == len(re.findall(r'\*', CIGAR_SUBSTITUTION_DELETION))