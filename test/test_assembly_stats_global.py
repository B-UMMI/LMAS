#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Module to test the LMAS script to get stats from an assembly file (FASTA)

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
ASSEMBLY_TEST =  "test/data/assembly.fasta"

def test_fasta_iter():
    fasta_iterator = utils.fasta_iter(ASSEMBLY_TEST)
    first_header, first_contig = fasta_iterator.__next__()
    
    assert first_header == 'k141_4'
    assert len(first_contig) == 1061
    assert sorted(set(first_contig)) == ['A', 'C', 'G', 'T']

def test_get_contig_lists():

    # test LMAS
    fasta_iterator = utils.fasta_iter(ASSEMBLY_TEST)
    contigs_len, contigs_len_over_1000, Ns_all, Ns_over_1000 = assembly_stats_global.get_contig_lists(fasta_iterator, 1000)

    # get quast result 
    with open('test/data/quast/report.tsv') as f:
        quast_report = list(csv.reader(f, delimiter='\t'))
        quast_contigs_len = int(quast_report[1][-1]) #second row, last col
        quast_total_len = int(quast_report[7][-1])
        quast_max_contig = int(quast_report[14][-1])
        quast_ns = float(quast_report[21][-1])
    
    assert len(contigs_len) == quast_contigs_len == 505
    assert sum(contigs_len) == quast_total_len
    assert max(contigs_len) == quast_max_contig

    assert len(contigs_len_over_1000) == 505
    assert min(contigs_len_over_1000) > 1000
    assert len(contigs_len) >= len(contigs_len_over_1000) 

    assert Ns_all >= Ns_over_1000
    assert Ns_all == quast_ns

def test_get_Nx():
    fasta_iterator = utils.fasta_iter(ASSEMBLY_TEST)
    contigs_len, contigs_len_over_1000, Ns_all, Ns_over_1000 = assembly_stats_global.get_contig_lists(fasta_iterator, 1000)
    n50_contigs = utils.get_Nx(contigs_len, 0.5)
    n50_contigs_over_min_len = utils.get_Nx(contigs_len_over_1000, 0.5)

    # get quast result 
    with open('test/data/quast/report.tsv') as f:
        quast_report = list(csv.reader(f, delimiter='\t'))
        quast_n50 = int(quast_report[17][-1]) #second row, last col

    assert n50_contigs == quast_n50 == 173935
    assert n50_contigs_over_min_len == 173935

    assert Ns_all == Ns_over_1000 == 0