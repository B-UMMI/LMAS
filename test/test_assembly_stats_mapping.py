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
import json
import pytest
from contextlib import contextmanager
from itertools import groupby
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
REFERENCE_TEST = "test/data/triple_reference.fasta"


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


def test_parse_paf_files():

    references = (x[1] for x in groupby(
        open(REFERENCE_TEST, "r"), lambda line: line[0] == ">"))
    alignment_dictitonary_list = assembly_stats_mapping.parse_paf_file(
        MAPPING_TEST, REFERENCE_TEST)

    assert len(alignment_dictitonary_list) == len([_ for _ in references])/2

    for alignment_dict in alignment_dictitonary_list:

        assert alignment_dict['Longest_Alignment'] <= alignment_dict['Reference_Length']

        sum_covered_bases = assembly_stats_mapping.get_covered_bases(
            alignment_dict['Covered_Bases'], alignment_dict['Reference_Length'])
        assert sum_covered_bases <= alignment_dict['Reference_Length']

        # TODO check with Mário
        #assert sum(alignment_dict['Alignment_Blocks']) <= alignment_dict['Reference_Length']
        #assert sum(alignment_dict['Alignment_Blocks']) >= sum_covered_bases

        for header in references:
            reference_name = header.__next__()[1:].strip()
            seq = "".join(s.strip() for s in references.__next__())
            if reference_name == alignment_dict['Reference']:
                assert alignment_dict['Reference_Length'] == seq/3


def test_get_mapping_stats():

    df = utils.parse_assemblies(
        "pytest_sample", "pytest_assembler", ASSEMBLY_TEST, MAPPING_TEST)
    
    alignment_dictitonary_list = assembly_stats_mapping.parse_paf_file(
        MAPPING_TEST, REFERENCE_TEST)

    to_plot_nax, to_plot_ngx, to_plot_lx, to_plot_phred, to_plot_coverage, json_dic = assembly_stats_mapping.mapping_stats(
        "pytest_sample", "pytest_assembler", df, alignment_dictitonary_list, 0.5, 0.9)

    # test JSON DICT for report
    assert list(json_dic.keys()) == ['sample_id', 'ReferenceTables']
    assert json_dic['sample_id'] == "pytest_sample"
    for reference in json_dic['ReferenceTables']:
        print(json_dic['ReferenceTables'][reference])
        assert json_dic['ReferenceTables'][reference]['assembler'] == "pytest_assembler"

        # test COMPASS
        assert 0 < round(json_dic['ReferenceTables'][reference]['breadth_of_coverage'], 4) <= 1
        #assert 0 < round(json_dic['ReferenceTables'][reference]['validity'], 4) <= 1 #TODO Still failing!
        #assert 0 < json_dic['ReferenceTables'][reference]['multiplicity']
        #assert 0 < json_dic['ReferenceTables'][reference]['parsimony'] # TODO <= 1? 

        # test Identity
        assert 0 < json_dic['ReferenceTables'][reference]['identity'] <= 1
        assert json_dic['ReferenceTables'][reference]['lowest_identity'] <= json_dic['ReferenceTables'][reference]['identity']

        # test contiguity
        assert json_dic['ReferenceTables'][reference]['contiguity'] >= 0
        if json_dic['ReferenceTables'][reference]['contiguity'] == 0:
            assert json_dic['ReferenceTables'][reference]['aligned_contigs'] == 0 
            assert json_dic['ReferenceTables'][reference]['aligned_basepairs'] == 0

        # conditional testing
        if json_dic['ReferenceTables'][reference]['breadth_of_coverage'] > 0:
           assert json_dic['ReferenceTables'][reference]['aligned_contigs'] > 0 
           assert json_dic['ReferenceTables'][reference]['aligned_basepairs'] > 0 
        elif json_dic['ReferenceTables'][reference]['breadth_of_coverage'] == 0:
            assert json_dic['ReferenceTables'][reference]['aligned_contigs'] == 0 
            assert json_dic['ReferenceTables'][reference]['aligned_basepairs'] == 0 


    # test nax df
    print(to_plot_nax.columns)