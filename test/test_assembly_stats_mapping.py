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
from templates import assembly_stats_mapping
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

    assert len(mapped_contigs.keys()) == 505  # number of mapped contigs
    ref_set = [item for sublist in list(
        mapped_contigs.values()) for item in sublist]
    assert sorted(set(ref_set)) == ['Bacillus_subtilis', 'Enterococcus_faecalis', 'Escherichia_coli', 'Escherichia_coli_plasmid', 'Lactobacillus_fermentum', 'Listeria_monocytogenes',
                                    'Pseudomonas_aeruginosa', 'Salmonella_enterica', 'Staphylococcus_aureus', 'Staphylococcus_aureus_plasmid1', 'Staphylococcus_aureus_plasmid2', 'Staphylococcus_aureus_plasmid3']


def test_parse_assemblies():

    df = utils.parse_assemblies(
        "pytest", "pytest", ASSEMBLY_TEST, MAPPING_TEST)

    assert df.shape == (516, 7)
    assert sorted(df.columns) == [
        '#N', 'Assembler', 'Contig', 'Contig Len', 'Mapped', 'Sample', 'index']
    assert len(df.Mapped.unique()) == 12

"""
def test_parse_paf_files():

    df = utils.parse_assemblies(
        "pytest", "pytest", ASSEMBLY_TEST, MAPPING_TEST)

    to_plot_nax, to_plot_ngx, to_plot_lx, to_plot_phred, json_dic = assembly_stats_mapping.parse_paf_files("pytest", df,
                                                                                                           MAPPING_TEST, 
                                                                                                           "test/data/triple_reference.fasta", "pytest",
                                                                                                           50, 50)
    print(json_dic)
"""