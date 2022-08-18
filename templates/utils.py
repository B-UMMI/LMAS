#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from itertools import groupby
import pandas as pd
import re
import logging

COLUMNS = ['Sample', 'Assembler', 'Contig', 'Contig Len',
           'Mapped', '#N']  # columns for dataframe

# colors for each assembler
COLOURS = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c',
           '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ebdb75', '#b15928']

ASSEMBLER_NAMES = ["ABySS", "GATBMiniaPipeline", "MetaHipMer2", "MINIA", "MEGAHIT", "metaSPAdes", "Unicycler", "SPAdes",
                   "SKESA", "VelvetOptimiser", "IDBA-UD"]

ASSEMBLER_PROCESS_LIST = ["ABYSS", "GATBMINIAPIPELINE", "MINIA", "METAHIPMER2", "MEGAHIT", "METASPADES", "UNICYCLER", "SPADES",
                          "SKESA", "VELVETOPTIMISER", "IDBA"]


def get_logger(filepath, level=logging.DEBUG):
    """

    :param filepath:
    :param level:
    :return:
    """
    # create logger
    logger = logging.getLogger(os.path.basename(filepath))
    logger.setLevel(level)
    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(level)
    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    # add formatter to ch
    ch.setFormatter(formatter)
    # add ch to logger
    logger.addHandler(ch)

    return logger


def log_error():
    """Nextflow specific function that logs an error upon unexpected failing
    """

    with open(".status", "w") as status_fh:
        status_fh.write("error")


def get_assember_name(assembly_file):
    """
    get assembler name from filename. Expected format: `XX_<AssemblerName>.fasta`
    :param assembly_file: path
    :return: filename
    """
    return os.path.basename(assembly_file).split('.')[0].rsplit('_')[-1]


def fasta_iter(fasta_name):
    """
    Given a fasta file. yield tuples of header, sequence.
    Modified from Brent Pedersen's Correct Way To Parse A Fasta File In Python.
    :param fasta_name: string with fasta file to parse
    :return: tuples with header, sequence (yield)
    """
    fh = open(fasta_name)

    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip().split(' ')[0]

        # join all sequence lines to one.
        try:
            seq = "".join(s.strip() for s in faiter.__next__())
        except StopIteration:
            print(headerStr)

        yield headerStr, seq


def get_mapped_contigs(paf_file):
    """
    Gets list with the sizes of the mapped contigs.
    In the paf file, the first col is the contig name,
    the second is the contig size (excludes gaps)
    :param paf_file: path to the PAF file
    :return: list with contig sizes
    """
    with open(paf_file) as f:
        mapped_contigs = [line.split()[0] for line in f]
    return mapped_contigs


def get_mapped_contigs_with_ref(paf_file):
    """
    Gets a dictionary with references and the the mapped contigs.
    In the paf file, the first col is the contig name,
    the sixth is the reference name
    :param paf_file: path to the PAF file
    :return: dict with reference and aligned contigs
    """
    mapped_contigs = {}
    with open(paf_file) as f:
        for line in f:
            if line.split()[0] not in mapped_contigs.keys():
                mapped_contigs[line.split()[0]] = [line.split()[5]]
            else:
                mapped_contigs[line.split()[0]].append(line.split()[5])
    return mapped_contigs


def parse_assemblies(sample_id, assembler, assembly, mapping):
    """
    Parses fastas and paf files and returns info on 'Assembler','Contig', 'Contig Len', 'Mapped' as dataframe
    :param sample_id, str with sample id
    :param assembler: str with assembler name
    :param assembly: assembly files
    :param mapping: paf files
    :return: pandas dataframe
    """
    df = pd.DataFrame(columns=COLUMNS)

    mapped_contigs = get_mapped_contigs_with_ref(mapping)

    fasta = fasta_iter(assembly)
    for header, seq in fasta:
        Ns = len(re.findall("N", seq.upper()))
        if header in mapped_contigs.keys():
            for reference in mapped_contigs[header]:
                is_mapped = reference
                df = df.append({'Sample': sample_id, 'Assembler': assembler, 'Contig': header, 'Contig Len': len(seq),
                        'Mapped': is_mapped, '#N': Ns}, ignore_index=True)
        else:
            is_mapped = 'Unmapped'
            df = df.append({'Sample': sample_id, 'Assembler': assembler, 'Contig': header, 'Contig Len': len(seq),
                            'Mapped': is_mapped, '#N': Ns}, ignore_index=True)

    df = df.reset_index()

    return df


def get_Nx(alignment_lengths, target):
    """
    Calculate NAx (x=target) form a list of contig lenghts
    :param alignment_lengths: list of aligned contig length sizes (unordered)
    :param target: percentage of total genome length, in float, from 0 to 1 (float)
    :return: na50 of the aligned contigs (also called NA50
    """
    sorted_lengths = sorted(
        alignment_lengths, reverse=True)  # from longest to shortest
    total_length = sum(sorted_lengths)
    target_length = total_length * target
    length_so_far = 0
    n50 = 0
    for contig_length in sorted_lengths:
        length_so_far += contig_length
        if length_so_far >= target_length:
            n50 = contig_length
            break
    return n50


def get_NGx(alignment_lengths, reference_length, target):
    """
    Calculate NGx (x=target) form a list of contig lenghts
    :param alignment_lengths: list of aligned contig length sizes (unordered)
    :param reference_length: genome gize
    :param target: percentage of known genome size, from 0 to 1 (float)
    :return: nx of the aligned contigs
    """
    sorted_lengths = sorted(
        alignment_lengths, reverse=True)  # from longest to shortest
    target_length = reference_length * target
    length_so_far = 0
    nx = 0
    for contig_length in sorted_lengths:
        length_so_far += contig_length
        if length_so_far >= target_length:
            nx = contig_length
            break
    return nx


def get_Lx(alignment_lengths, ref_len, target):
    """
    Returns the number of contigs, ordered by length, that cover at least 'target'% of the reference sequence.
    :param alignment_lengths: list with length of mapped contigs for the reference
    :param ref_len: int with the expected reference length
    :param target: target % of the reference sequence for Lx metric, from 0 to 1 (float)
    :return: int with the number of contigs that represent
    """
    sorted_lengths = sorted(
        alignment_lengths, reverse=True)  # from longest to shortest
    target_length = ref_len * target

    if sum(sorted_lengths) < target_length:
        return 0
    
    if len(sorted_lengths) == 1:
        if sorted_lengths[0] >= target_length:
            return 1

    length_so_far = 0
    Lx = 0
    for contig_length in sorted_lengths:
        length_so_far += contig_length
        if length_so_far <= target_length:
            Lx += 1
    return Lx


def is_number(n):
    """
    Verify if n is a number by trying to set it to float
    :param n: value to be verified
    :return: Bool is value is a number or not
    """
    try:
        float(n)
        return True
    except ValueError:
        return False


def parse_cs(string):
    """
    Parses PAF's cigar string to obtain the number of snps and indels.
    The sybstitutions are marked with "*" followed by the reference base and the substituting base,
    the insertions with '+' followed by the bases inserted in comparison to the reference, and
    the deletions with '-' followed by the bases deleted in comparison to the reference.
    :param string: Cigar-like string to be parsed
    :returns
    """
    indel = []

    exact_matches = sum(map(int, re.findall(r':([\d]+)', string)))

    # substitutions
    snps = len(re.findall(r'\*', string))

    # insertions
    insertions = re.findall(r'\+([a-z]+)', string)
    for insertion in insertions:
        indel.append('+' + str(len(insertion)))

    # deletion
    deletions = re.findall(r'-([a-z]+)', string)
    for deletion in deletions:
        indel.append('-' + str(len(deletion)))

    return exact_matches, snps, indel


def cs_parse_deletion(string):
    """
    Parses PAF's cigar string to obtain the length of the deletion
    """
    deletion_lenght = 0

    deletions = re.findall(r'-([a-z]+)', string)
    for deletion in deletions:
        deletion_lenght += len(deletion)
    
    return deletion_lenght

def cs_get_matched_bases(string):
    return sum(map(int, re.findall(r':([\d]+)', string)))


def adjust_reference_coord(coord, ref_len):
    """

    :param coord:
    :param ref_len:
    :return: int coord adjusted to the real length of the reference
    """
    if coord <= ref_len:
        return (coord)
    elif coord <= 2 * ref_len:
        return (coord - ref_len)
    else:
        return (coord - (2 * ref_len))


def check_overlap(list_of_coords):
    """
    Function that takes a list of coords and checks if there is an overlap
    :param list_of_coords: list of lists containing start (0 based, closed) and stop (0 based, open) positions
    return: True if there is an overlap, False otherwise
    """
    sorted_list_of_coords = sorted(list_of_coords, key=lambda x: x[0])
    list_of_ranges = []
    for coods in sorted_list_of_coords:
        x = range(coods[0],coods[1])
        list_of_ranges.append(set(x))
    for i in range(len(list_of_ranges)-1):
        if len(list_of_ranges[i].intersection(list_of_ranges[i+1])) > 0:
            return True
    return False

def get_check_overlap(list_of_coords):
    """
    Function that takes a list of coords and returns a list with overlap lengths 
    :param list_of_coords: list of lists containing start (0 based, closed) and stop (0 based, open) positions
    return: list with overlap lengths 
    """
    sorted_list_of_coords = sorted(list_of_coords, key=lambda x: x[0])
    list_of_ranges = []
    overlaps = []
    for coods in sorted_list_of_coords:
        x = range(coods[0],coods[1])
        list_of_ranges.append(set(x))
    for i in range(len(list_of_ranges)-1):
        overlaps.append(len(list_of_ranges[i].intersection(list_of_ranges[i+1])))
    return overlaps