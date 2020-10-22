#!/usr/bin/env python3
import os
import fnmatch
from itertools import groupby
import pandas as pd
import re
import logging

COLUMNS = ['Assembler', 'Contig', 'Contig Len', 'Mapped']  # columns for dataframe

# Dic for pretty print of reference names
REFERENCE_DIC = {
    "BS.pilon.polished.v3.ST170922": "Bacillus subtilis",
    "Enterococcus_faecalis_complete_genome": "Enterococcus faecalis",
    "Escherichia_coli_chromosome": "Escherichia coli",
    "Lactobacillus_fermentum_complete_genome": "Lactobacillus fermentum",
    "Listeria_monocytogenes_complete_genome": "Listeria monocytogenes",
    "Pseudomonas_aeruginosa_complete_genome": "Pseudomonas aeruginosa",
    "Salmonella_enterica_complete_genome": "Salmonella enterica",
    "Staphylococcus_aureus_triple_chromosome": "Staphylococcus aureus"
}


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
        headerStr = header.__next__()[1:].strip().split()[0]

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
    Gets list with the sizes of the mapped contigs.
    In the paf file, the first col is the contig name,
    the second is the contig size (excludes gaps)
    :param paf_file: path to the PAF file
    :return: dict with contig names and matching reference
    """
    mapped_contigs = {}
    with open(paf_file) as f:
        for line in f:
            mapped_contigs[line.split()[0]] = line.split()[5]
    return mapped_contigs


def parse_assemblies(assemblies, mappings):
    """
    Parses fastas and paf files and returns info on 'Assembler','Contig', 'Contig Len', 'Mapped' as dataframe
    :param assemblies: list of assembly files
    :param mappings: list of paf files
    :return: pandas dataframe
    """
    df = pd.DataFrame(columns=COLUMNS)

    for fasta_file in assemblies:

        filename = get_assember_name(fasta_file)
        mapped_contigs = get_mapped_contigs(fnmatch.filter(mappings, '*_' + filename + '.*')[0])

        fasta = fasta_iter(fasta_file)
        for header, seq in fasta:
            if header in mapped_contigs:
                is_mapped = 'Mapped'
            else:
                is_mapped = 'Unmapped'

            df = df.append({'Assembler': filename, 'Contig': header, 'Contig Len': len(seq), 'Mapped': is_mapped},
                           ignore_index=True)

    df = df.reset_index()

    return df


def get_N50(alignment_lengths):
    """
    Callculate n50 form a list of contig lenghts
    :param alignment_lengths: list of aligned contig length sizes (unordered)
    :return: n50 of the aligned contigs (also called NA50
    """
    sorted_lengths = sorted(alignment_lengths, reverse=True)  # from longest to shortest
    total_length = sum(sorted_lengths)
    target_length = total_length * 0.5
    length_so_far = 0
    n50 = 0
    for contig_length in sorted_lengths:
        length_so_far += contig_length
        if length_so_far >= target_length:
            n50 = contig_length
            break
    return n50


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


def adjust_reference_coord(coord, ref_len):
    """

    :param coord:
    :param ref_len:
    :return: int coord adjusted to the real length of the reference
    """
    if coord <= ref_len:
        return coord
    elif coord <= 2 * ref_len:
        return coord - ref_len
    else:
        return coord - (2 * ref_len)
