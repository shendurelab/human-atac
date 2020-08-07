#!/usr/bin/env python
import argparse
import subprocess
import sys
import os
import gzip
import io
import itertools
import time
import json
import collections
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def get_row_col_count(index_set):
    if len(index_set) == 96:
        ncol = 12
        nrow = 8
    elif len(index_set) == 384:
        ncol = 24
        nrow = 16
    else:
        raise ValueError('%s indices in index set, but only 96 or 384 are supported when converting to well IDs.')
    
    return nrow,ncol

def chunk_list(l, n):
    """
    Chunks list into N sized chunks as list of list.
    """
    if n <= 0:
        raise ValueError('Chunk size of %s specified, which is invalid, must be positive int.' % n)

    results = []
    for i in range(0, len(l), n):
        results.append(l[i:i + n])
    return(results)

def merge_dicts(x, y):
    """
    Merges the entries specified by two dicts.
    """
    z = x.copy()
    z.update(y)
    return z

def pad_well_col(well_col, zero_pad, id_length):
    if zero_pad:
        template = '%%0%sd' % id_length
    else:
        template = '%s'
    col_id = template % (well_col)
    return col_id
    
def get_well_id(i, row_ordered, nrow, ncol, zero_pad_col, id_length):
    if row_ordered:
        well_row = chr(65 + int(i / ncol))
        well_col = (i % ncol) + 1
    else:
        well_row = chr(65 + (i % nrow))
        well_col = int(i / nrow) + 1

    well_id = '%s%s' % (well_row, pad_well_col(well_col, zero_pad_col, id_length))
    return well_id


def get_well_dict(index_set, row_ordered=True, zero_pad_col=True):
    """
    Transforms an set of indices into a well dict.
    
    Args:
        index_set: list of all indices
        row_ordered (bool): true if these indices are row ordered and false if column ordered
        zero_pad_col (bool):  true if you want A01 vs. A1 for example
        
    Returns:
        dict: dict elements in the index set to a well.
    """
    nrow, ncol = get_row_col_count(index_set)
    
    id_length = len(str(len(index_set)))
    mapping = dict()
    
    for i,element in enumerate(index_set):
        well_id = get_well_id(i, row_ordered, nrow, ncol, zero_pad_col, id_length)
        mapping[element] = well_id

    return mapping


def get_row_col_matrix(row_index_set, col_index_set):
    """
    Returns the resulting matrix when using a row and column of indices (list of tuples)
    row_index_set: the index set that would be taken from wells A-H (or equivalent) on that plate. Should be dim 8 or 16.
    col_index_set: the index set that would be taken from wells 1-12 (or equivalent) on that plate. Should be dims 12 or 24.
    
    """
    if not ((len(row_index_set) == 8 and len(col_index_set) == 12) or (len(row_index_set) == 16 and len(col_index_set) == 24)):
        raise ValueError('Unexpected row or column length. get_row_col_matrix expects the row index set to be the index set from A-H or equivalent (8 or 16 in length) and the col index set to be from wells 1-12 or equivalent (12 or 24 in length). Make sure that you have not swapped the two.')

    index_set = []
    for row_index in row_index_set:
        for col_index in col_index_set:
            index_set.append((row_index, col_index))
    return index_set

def get_pcr_plate_dict(p5_pcr_plate, p7_pcr_plate, zero_pad_col=True):
    """
    Returns a mapping for (p7_index, p5_index) to a well ID that maintains a constant well column for row/col combinations
    along with an indication of the original PCR plate row/col wells that were combined.
    """
    nrow, ncol = get_row_col_count(p7_pcr_plate)
    nrow_p5, ncol_p5 = get_row_col_count(p5_pcr_plate)
    if nrow != nrow_p5 or ncol != ncol_p5:
        raise ValueError('P7 and P5 plate must have equal dimensions: (%s, %s), (%s, %s)' % (nrow, ncol, nrow_p5, ncol_p5))

    p7_rows = chunk_list(p7_pcr_plate, ncol)
    p5_cols = chunk_list(p5_pcr_plate, nrow)

    # Build a map of each individual plate too to get the coordinates
    p5_plate_map = get_well_dict(p5_pcr_plate, row_ordered=False)
    p7_plate_map = get_well_dict(p7_pcr_plate, row_ordered=True)

    full_layout_mapping = dict()
    
    for row in range(0, len(p7_rows)):
        for col in range(0, len(p5_cols)):
            row_col_layout = get_row_col_matrix(p5_cols[col], p7_rows[row])
            row_col_layout_mapping = get_well_dict(row_col_layout, row_ordered=True, zero_pad_col=zero_pad_col)
            
            row_col_layout_mapping_final = {}
            for k,v in row_col_layout_mapping.items():
                k_p5, k_p7 = k
                row_col_string = '_row%s_col%s' % (p7_plate_map[k_p7], p5_plate_map[k_p5])
                row_col_layout_mapping_final[k] = '%s%s' % (v, row_col_string)

            full_layout_mapping = merge_dicts(full_layout_mapping, row_col_layout_mapping_final)
        
    return full_layout_mapping

def get_two_level_barcode_string(tagmentation_i7_seq, pcr_i7_seq, pcr_i5_seq, tagmentation_i5_seq, tagmentation_to_well, pcr_to_well, well_ids):
    """
    Helper function to construct final barcode string for two-level indexed Tn5 protocol.
    """
    if well_ids:
        tagmentation_well_string = tagmentation_to_well[(tagmentation_i5_seq, tagmentation_i7_seq)]
        pcr_well_string = pcr_to_well[(pcr_i5_seq, pcr_i7_seq)]
        outputs = [tagmentation_well_string, pcr_well_string]
        barcodes_string = '-'.join(outputs)
    else:
        outputs = [tagmentation_i7_seq, pcr_i7_seq, pcr_i5_seq, tagmentation_i5_seq]
        barcodes_string = ''.join(outputs)

    return barcodes_string

def get_barcode_string(tagmentation_i7_seq, pcr_i7_seq, pcr_i5_seq, tagmentation_i5_seq, tagmentation_i7_to_well, tagmentation_i5_to_well, pcr_to_well, well_ids):
    """
    Helper function to construct final barcode string for 3LV (or two-level ligation based) protocol.
    """
    if well_ids:
        tag_n7_string = tagmentation_i7_to_well[tagmentation_i7_seq]
        pcr_well_string = pcr_to_well[(pcr_i5_seq, pcr_i7_seq)]
        tag_n5_string = tagmentation_i5_to_well[tagmentation_i5_seq]
        outputs = [tag_n7_string, tag_n5_string, pcr_well_string]
        barcodes_string = '-'.join(outputs)
    else:
        outputs = [tagmentation_i7_seq, pcr_i7_seq, pcr_i5_seq, tagmentation_i5_seq]
        barcodes_string = ''.join(outputs)

    return barcodes_string

def correct_barcode(barcode, mismatch_map):
    """
    Correct an observed raw barcode to one of a list of whitelists of mismatches.
    Args:
            barcode (string): barcode sequence to be corrected
            mismatch_map (list of dict dict): list of dict of mismatched sequences to real sequences
    Returns:
            string: corrected barcodes or None if barcode not correctable.
    """
    for mismatch_whitelist in mismatch_map:
        corrected = mismatch_whitelist.get(barcode, None)

        if corrected:
            return corrected

    return None


def generate_mismatches(sequence, num_mismatches, allow_n=True):
    """
    Generate a list of mimatched sequences to a given sequence. Must only contain ATGC.
    This is heavily based on a biostars answer.
    Args:
        sequence (str): The sequence must contain only A, T, G, and C
        num_mismatches (int): number of mismatches to generate sequences for
        allow_n (bool): True to allow N bases and False if not
    Yield:
    """
    letters = 'ACGT'

    if allow_n:
        letters += 'N'

    sequence = sequence.upper()
    mismatches = []

    for locs in itertools.combinations(range(len(sequence)), num_mismatches):
        sequence_list = [[char] for char in sequence]
        for loc in locs:
            orig_char = sequence[loc]
            sequence_list[loc] = [l for l in letters if l != orig_char]

        for poss in itertools.product(*sequence_list):
            mismatches.append(''.join(poss))

    return mismatches


def construct_mismatch_to_whitelist_map(whitelist, edit_distance, allow_n=True):
    """
    Constructs a precomputed set of all mimatches within a specified edit distance and the barcode whitelist.
    Args:
        whitelist (set of str): set of whitelist sequences
        edit_distance (int): max edit distance to consider
        allow_n (bool): True to allow N bases and False if not
    Returns:
        dict: mapping of mismatched sequences to their whitelist sequences
    """

    mismatch_to_whitelist_map = [None] * (edit_distance + 1)

    mismatch_to_whitelist_map[0] = {k: k for k in whitelist}

    conflicting_mismatches = []  # tracks conflicts where mismatches map to different sequences

    # Doesn't really matter as correction function will never see it,
    # but exclude any perfect matches to actual seqs by mismatches
    conflicting_mismatches.extend(list(whitelist))

    for mismatch_count in range(1, edit_distance + 1):
        mismatch_to_whitelist_map[mismatch_count] = {}

        for sequence in whitelist:
            sequence = sequence.upper()

            # Generate all possible mismatches in range
            mismatches = generate_mismatches(sequence, num_mismatches=mismatch_count, allow_n=allow_n)

            # Construct a mapping to the intended sequences
            for mismatch in mismatches:
                # Check for conflict with existing sequence and track if so
                if mismatch in mismatch_to_whitelist_map[mismatch_count]:
                    conflicting_mismatches.append(mismatch)
                mismatch_to_whitelist_map[mismatch_count][mismatch] = sequence

        # Go back and remove any conflicting mismatches
        for mismatch in set(conflicting_mismatches):
            if mismatch in mismatch_to_whitelist_map[mismatch_count]:
                del mismatch_to_whitelist_map[mismatch_count][mismatch]

    return mismatch_to_whitelist_map


def reverse_complement(x):
    complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    xrev = x[::-1]
    xrevcomp = ''.join([complements[z] for z in xrev])
    return xrevcomp

def get_barcode_seqs(r1_name, nextseq, two_level_indexed_tn5):
    """
    Extract the correct sequences from the R1 name.
    """
    # In 3LV runs, this is the the P5 + P7 index seq with a + in between
    # which is 20 + 1 + 20 (so have to skip "+" position)
    # Similar for two-level indexed Tn5, but Tn5 barcodes are 8bp
    if not two_level_indexed_tn5:
        barcodes = r1_name[-41:]

        tagmentation_i7_seq = barcodes[0:10]
        pcr_i7_seq = barcodes[10:20]

        if nextseq:
            pcr_i5_seq = barcodes[31:41]
            tagmentation_i5_seq = barcodes[21:31]
        else:
            pcr_i5_seq = barcodes[21:31]
            tagmentation_i5_seq = barcodes[31:41]
    else:
        barcodes = r1_name[-37:]

        tagmentation_i7_seq = barcodes[0:8]
        pcr_i7_seq = barcodes[8:18]

        if nextseq:
            pcr_i5_seq = barcodes[27:37]
            tagmentation_i5_seq = barcodes[19:27]
        else:
            pcr_i5_seq = barcodes[19:29]
            tagmentation_i5_seq = barcodes[29:37]

    return tagmentation_i7_seq, pcr_i7_seq, pcr_i5_seq, tagmentation_i5_seq

def indexsplitter(indexrange):
	if len(indexrange) < 3:
		indexout = [int(indexrange)-1]
	elif "-" in indexrange or ',' in indexrange:
		range_list = [x for x in indexrange.split(",")]
		indexout = []
		for myrange in range_list:
			index_range = myrange.split('-')
			
			if len(index_range) == 1:
				start = int(index_range[0]) - 1
				end = start + 1
			elif len(index_range) == 2:
				start = int(index_range[0]) - 1
				end = int(index_range[1])
			else:
				raise ValueError('Invalid index range %s' % myrange)

			indexout.extend(range(start, end))
	else:
		raise ValueError('Invalid format for index range: %s' % indexrange)
	return indexout

def required_index(a):
    """
    Helper function to take a list of index lists and return whether it needs to be included as an index in demultiplexing.
    """
    return len(set(tuple(a_i) for a_i in a)) != 1
    
def get_sample_lookup(samplesheet, pcri7, tagi7, tagi5, pcri5):
    """
    Takes a samplesheet file handle and returns a list of True/False indicating whether index was used in lookup and a lookup of tuples of index combinations to sample names.
    Args:
        samplesheet (file handle): opened file handle to samplesheet
        tagi5 (list): list of tag/lig i5 indices
        pcri5 (list): list of pcr i5 indices
        pcri7 (list): list of pcr i7 indices
        tagi7 (list): list of tag/lig i7 indices
    Returns:
        (list of bool, dict of tuple to sample name): [use_pcri7, use_tagi7, use_tagi5, use_pcri5] where each entry indicates usage of that barcode in indexing and lookup of (tagi7,tagi5), for example.
    """
    pcri7_indices = []
    tagi7_indices = []
    tagi5_indices = []
    pcri5_indices = []
    samples = []
    for line in samplesheet:
        if line.startswith('sample_id\tranges'):
            continue

        entries = line.strip().split()
        sample, indices = entries
        indexsplit = indices.split(':')

        samples.append(sample)
        pcri7_indices.append(indexsplitter(indexsplit[1]))
        tagi7_indices.append(indexsplitter(indexsplit[0]))
        tagi5_indices.append(indexsplitter(indexsplit[3]))
        pcri5_indices.append(indexsplitter(indexsplit[2]))

    use_pcri7 = required_index(pcri7_indices)
    use_tagi7 = required_index(tagi7_indices)
    use_tagi5 = required_index(tagi5_indices)
    use_pcri5 = required_index(pcri5_indices)

    index_mask = [use_pcri7, use_tagi7, use_tagi5, use_pcri5]
    index_whitelists = [pcri7, tagi7, tagi5, pcri5]
    index_lists = [pcri7_indices, tagi7_indices, tagi5_indices, pcri5_indices]
    
    sample_lookup_table = {}
    for sample_index,sample in enumerate(samples):
        indices_to_use = []
        for i,(use_index,index_list) in enumerate(zip(index_mask, index_lists)):
            if use_index:
                indices_to_use.append([index_whitelists[i][index_i] for index_i in index_list[sample_index]])

        for combination in list(itertools.product(*indices_to_use)):
            sample_lookup_table[combination] = sample
        
    return index_mask, sample_lookup_table

# Set up constants
## 96-well 3LevelV2 and well mappings
lig_i7_list = ["TTAGTATATA","TGCCTCGCAA","GACTACGGCC","CGGCCGCTCC","TTCATCGTAA","GAGATTCCAA","CGCCTACCAA","GTTATATGCA","TTACGAATGG","ATTCTGCGGT","GCTGATAATT","AACGCAATGG","CCAACCGCGA","CTTCTCGACT","AATGATGCTC","CGGTCCGCGG","GCAAGTTAGC","AACCAGCATG","CCTCTAGGAC","AGAGGACTGG","GCTCTTACTT","CCTTGGCTCT","TGGCGTTGGC","TACTGGACGG","CCTTAACGCC","GTCGAGGTAA","GATGATTGGA","GAACTTGCGG","TAGATCCTAA","GCCGACGATA","AGTTACGCCG","GAGGAGAGCT","TACTCGAACT","TTATAGGAGA","CGGAGCGTCA","GCGACTCGAT","TATCGGCAGT","GAGTATACCT","TGGAGGCGAA","CGGATAAGCT","GGCTGCGACG","CTGAAGTCAA","ACTAGACGTT","TTGGTCGAAT","CTAACCTACC","CTGCGGACGT","CTGCAGATCC","CTTCTCTTAG","ATTCAATTAA","TCTGCCGGCA","GATTCGGTCA","ACCAGCGGTT","ATGACGGCTA","AAGCTATGCC","GTAAGGCCGA","TTCTTGACCG","TCTAGCAAGT","AGGCAGTACC","AATTCAGAAC","TCTGGACTCA","CAATTCGTTA","ATGGCAACTT","CTCGACCGGC","TACTCCTCGA","GACTAGCATT","CTCTGGAGAC","AAGTTCCAGC","ACCTTCTGGA","AAGGATCCAA","TACCGCTCCT","TCTGCTATAA","CGAGCGGCCT","TGCTGCATCC","GCTCCAACGC","TAATGGTAGA","TGGCTGACCA","GTTCTTATAA","GCTATCCTCC","ACGATAGCTG","GATCCGCCAT","TAGATAATGG","CCTATTACCT","AGAAGCGATC","ACTACCGCGC","CGAGGCCGGA","TATTGCAGCT","AGAATTAGCA","TGATCCTGGC","GCCGTAATCG","AAGGCGTCGC","TGCATGGCCA","CAGCCATGAT","CGTAAGCGAA","GTCGGCAAGA","GCAACTCCGG","ATATAGGTTG"]
pcr_i7_list = ["GAATGCCTTG","CGGACTGGCC","CTGGCAGCGG","TCATCAGCCA","TAACCAGTTA","AGCTCATTCG","CAACGGTCTT","TGGTATCAGA","GGTAACTACG","CTTCAGGCCA","ATCCTTCAAC","TGATTGAATG","TGCTTAGCGA","CGAGTTCGCC","CGGCATCTTG","TCAACGAGGT","GGTCTAGGAA","AAGTCAACGG","TAACCAATCG","ATCGGCTTGA","CAATCTCAAT","CAGTTGCTTG","ATGCTAACCT","CGGATCTCCG","TGACGCGACC","TCTGACGAAC","GATCATGATA","TTCTCTACGA","TACGCAGGTT","TGCGAATCGG","CTAAGAGTTA","TCCATACCTG","ACCGGTCTGC","TAGAACCTCC","TGACCATCAG","GACTCGAGGA","GCATTAGGCG","GGACTTATGG","TTGGTTCTCA","TAATTCCGGT","CTACTGCAAG","CGTAGCCGTA","TTGACTTCGT","TCAGCGCAAT","ATATGGTACC","TGCCGGTACC","GGTACTTCCA","AATTGAGTTA","GCGGCGCAAG","AATACCTGCC","AGCATACGGC","GCCTCGGTAA","ACGGACCTGC","CATTCTGATG","GCTTACTATA","TCCAATGGTT","ACGGTAATTG","TCGGAGTATG","GGTTAGTTGG","CTGCATTACG","TGCGGCCTGG","CATCTATTCG","GCAGCTGATT","ATATCTTCCG","GTCGTACGCG","TTGATGCCTC","TTGGAACTGG","GGAGAATGGC","TTCCATTCTT","CCGAAGGTTA","CCGGCTGAGC","ACGAGGAGCC","AAGATCAGTT","ATCCGAACGG","CTGCGCTGGT","TGCCATGGAA","AGCTGAACGC","GAGACGTTCT","CTGACTACTT","CGCCGATATC","TGGATCGTTA","GTACCTGAAT","GAATTCTCCA","AGACCAGGTT","CTAGCTTCTT","GGCAACCTTC","CCGTTAGCAA","GGCGTCTGCC","CGCGCAACCG","AGTCTTAATT","TAAGGAACGG","TCGTTCATTA","ACTCCAGATT","CGATATATCT","GGCTGGCTCT","GGCGGATACC"]
pcr_i5_list = ["GACCTCCTTG","ACCAGTTCAG","TTCTGGCGCA","GGCAATCAAG","TTGGCCAGGT","TCTCCTCCTG","CTTCAGTAGT","AAGGTTAAGT","TGAGCCGCGG","TGGCCGTCTT","AACCTTCTCT","GGCGCCATCG","GTCGACCATT","CATTACCAAG","CTAGTTAAGT","CATATTGCAT","GCCAAGGCAA","GATCCATAAC","ATGAATACCA","TTATTGCTAA","TATTCCTTGC","GGTACCGGCT","GTACTGGTTG","CGTTGATGAT","ATGCCTCTTC","TTATATTGAA","GTTCGTCTGA","ATCCGGATAA","TTGCCTTATC","ATCTAGTTAG","AGATAGAGGT","TTCTAATCCG","AAGATGGATT","TTCCAAGTAA","ATGAAGCTTG","ATCCTTGGAG","GTTCCGTATC","GCCGGAGCGG","AGATGGCGCA","TTATTAAGAA","CTCTTCCAGG","TTAGCGGTTG","AGTAGATCTA","GTATCAATAT","GGTCTTCGGT","CCTAGTCTAT","GACCGTTACT","TATGCTTACT","CAGGAATCGG","AGGAGTAAGG","TGAGGTTAGC","GCATCGCGCA","GAACGTCAGG","TGAGTTGGAC","ATAGCATCAA","TCTAGTTCCA","TGCGTCTATA","ATAGAGAGCC","TCTTAATCAG","CGGAAGACTT","GCCAGACTCA","TAACTGGCTT","AAGGCCATCA","GTATACGTAA","CCAATTCCAT","AATATAGCGG","GTTCGAGGCG","ACGTACTAGG","ATCGAATCTG","TCTGCGTCGC","TAAGCTCGCC","GGCATAACCG","AGTCGCGTCG","TGGACCAAGG","CTGATCAGAG","TCTCATTGCC","ACTCGGCATA","GAGCTCAGCC","ATTGAGGATA","TTAGCGGAAC","CCTTCATCCG","GGCAGCAGTT","TAAGTAAGTA","TATCAAGGTC","GCGGTTGGAA","GAATGAATAA","GAAGAGTATT","TACCTTAGCT","CGTCAGTCAA","GCGCCTAGTT","GTAGTAGTCC","ACGCTTCGTC","GATAGCCGAT","AGTCATGGAT","TCTTAAGGAC","CCGGTCCTAA"]
lig_i5_list = ["TGCTCCTTAT","CTCGATCTTC","GCAGGTCGTC","ACTAATCCAG","GCTGGTTATG","ATGAGTAATA","CAATAATTAA","GAAGCAGCGG","TGCGCCGAAG","TCCTGGTCAA","GTTAATTCTA","GCGTAAGATC","CTCAGAATAA","AACCTATAGT","ACCGGAAGAA","ATTGCATAGT","ATTCCTACCG","ACGCGACGGC","GTTCTCTCCT","TTATAGCTCT","AATCGTTGAT","TATTACTCTA","ATGCCAGCAA","CAAGTAGGAC","CTTACGTTGA","TTCGCGATGG","TCGATAAGAA","TTAGAGCTGA","TAGGTCCGCA","ATGCGATTGA","TACTATTAGA","CTATGGCCTA","CCGAACTATG","ACGATCAAGG","GGATGGTATT","GGAGGCCTCC","TTAACGACCG","AGACCGCAAG","ATTATCGATT","AGGTTGCAAG","TGCCGCAGAG","TTACTAGTAC","ACTTGCTGAT","TGAACGGCGT","TAAGGCTGGT","AGGCCGGCCA","AATCTCGCGT","CCGCCTCCGT","CGCTACGATT","GCTTGAGTCT","GAACTAACCT","GCGTCCAGGT","GGTCGACTTG","ACCATAATAA","GTCCTATGAA","GCTACCGTAG","TGGCTCATAT","TATCTGCCTA","CATGCGAGAA","ACTGGTAGGA","CCATATGCTC","GGATTCTCGG","TAGTCGGATG","TGGTCTTGAA","CGTCGCCTGC","GTTCGGAGAT","TCAGATTCTT","TATATGACGT","CCGTCGAATG","CAGCGTCCAA","TCGTTGCGGC","AAGAAGTCTA","ATAATCGCGG","ATGGAGCTAC","GCAGAGCCAT","ACTCTCTAGT","CAATGAGAAG","CTAATATGAT","TCGATGATAC","TTAACTCATA","CTACGACTAT","CCGGCGTTCG","AATATATCAA","AGACTACTAG","GGTTGGCGGT","GAATGAGGTT","TGCAGAGATG","CGTATTATAC","GTTAGGAGCG","CAAGCATACC","CCTATTGAGT","TGGTATTATT","GAATTATTCG","AAGCGGCATT","CTTGCGTAAC","CCTGCGTATT"]

lig_i7_to_well = get_well_dict(lig_i7_list, row_ordered=True)
pcr_to_well = get_pcr_plate_dict(pcr_i5_list, pcr_i7_list)
lig_i5_to_well = get_well_dict(lig_i5_list, row_ordered=False)

lig_i7 = set(lig_i7_list)
pcr_i7 = set(pcr_i7_list)
pcr_i5 = set(pcr_i5_list)
lig_i5 = set(lig_i5_list)

## 384-well 3LevelV2 and well mappings
lig_i7_list_384 = ["TTAGGCGTCC","GGTCTAACTG","ATTGAGACGG","AACGTACCAT","GTTCAGTTAC","GTTAACCTTA","GGAGTAGTAG","CAGGTTGAGA","AACCAATAAC","GGTACCTATT","AGCTGCCTCA","TTGGTTGGTA","CTGATATCGG","TTGATGGAAC","GGACCAGAAG","ACCTGAGGTC","TTGCAATGAG","ACTAACTGAC","GGTAAGAGGT","GTAATCGCAA","ATATGATGAA","CCTACGGAAG","GACTCTCCGA","ATAGATACGT","TATGAAGCAA","GTAGAGACGA","TTCGCATAAC","TAGATTCGCC","TTCTTCGCGG","TTCTTGGTCT","GTATAGATTA","GACTATGACT","GCTCCGCGAA","GGAACGATTC","GTACCTTCGT","TTATGCGACT","AACGGTTGGT","CTAGAAGGAA","GTATGGAGGA","CGTAAGTAGC","TAGTCGTCAA","GTAGTATGGA","AATGGACGTA","CCGTTCGCTG","CAGTCAATAT","GTCGGACTGA","TTCCATGCGC","GGACGGCATG","GGAGCAACGT","CGTCCTAGCT","ATATTAGTAG","AGTTCCTTCT","GCCTTCAAGG","TAACGCCTCA","GAATACGTCT","GATGAGCCTT","TACGGAGACT","AGAGGTCATT","CGGCCAGTTA","TAGCGCTTAA","CAACCTTATG","CTTCGACGAA","AGAGTCCGAG","CTATGGTTCG","AACTATTGGC","CTTGCGAGGT","TTCGGCTCGT","CTCAAGTCTG","ATTAGGCCTG","TGGCTAACGT","AGCTAATATT","AGGAATCTTC","AAGACCAGAC","ACCATAAGCG","TGCTTATTAC","ATAGAGCATT","AATTATCTTG","GTTCCGCGCG","CCGAAGTCCG","TTACCATGAT","TATAGATGGC","GCATATCTGG","ACTCTAGAAG","CAAGAATATG","CTGCTGCGCG","CATATTACGG","CGGCGGAATG","TGCTTCTGAA","GTATTCTTCT","CAGTAGTAAT","CCGCGGAGAC","GCTCGTCGGC","CTTAGAGGCA","GACCAAGACG","AGCAGATGGT","AATTAACCAT","CAGACTCGCT","CTTCGCCGTT","CCTGCCAACC","GAGATCGCGG","AATTGCCGTC","GCGCGGTCCG","GCGAAGGCTC","CATCGACTCA","CTCGGTTGCA","GTTAAGTATT","AGTTCCGGAC","CAGCGGTCGG","GACGCGGTAG","TTCTTCTATT","CAATTATAGT","ATAGGAGTTC","GTTATTAATT","ATGATCGCCA","ATCATCAGAA","GTTCGGCCTT","TTACGTCTAG","GCTGCATAGC","TGATTGCTTC","TAACTTGGAG","AACTTACGCT","GGTTCGTACC","CGCTCCTTGG","GCTTGCGCGT","AATCTAGTCC","TTGGAGAACG","GGTAGTCTTA","TGGAACCGTA","GAGTTAAGGT","GTTATCGGTA","CTACGAGTCG","CCTAATTCTG","TGCCTTGACG","TCTAATTATT","CGCAACCTGG","AGGCGGTCTC","TCTCGTTGAG","AATTACTAAT","CTATCCGGAG","TTCGTTCCTG","CGACCTATAC","CTGGCTTAGT","TACGCCGTAT","GAGTCCTGCA","TGGCGCAACT","CCGTCTACTA","AGGCTGATTA","CTCCATAATC","ATTACGTAGG","ACGTATGCGC","GTCTAATAGG","TCAGAACGAC","GTCCAGTCGT","TGATGCTGAT","CGCGTTACCT","TCTCTGATCG","GGTATCATCT","GCCTGGCAGT","TGGTTACGCT","CCTACCGTCC","CGCTCCAGAA","TTCCGCGTTA","CCGACGCGCC","GATCTTATCC","CGCGAGGATG","TGGAATATAG","GCGGCTACTC","GTTATGAGAA","CAACGTTGCC","AGTAGCCTTG","GCCGGTTAAC","AAGATATAGA","TCTCAGCGGC","TCCAAGGAGA","GACCAATTCT","GGACCTTGCA","TGATCAATTA","TTCTCCGCTA","TTGGACGCGG","AACTGCTAAG","CGAGTCGAAG","TCAGGATTAA","AAGGACCTTA","TCCTCTTAAC","TCATAGGCTT","GTTACGACTC","GAGGTCCAGT","GGAAGTCTGG","CATGATCGTT","CTCATGATTG","GGCAGAGGAG","GTTCTAGCCA","CAGCGGATAA","GCGAGAGAAG","TGAATAACTA","TGGAGAGGCC","GAGGACCGAT","CGCTGGCAAG","CTTGATTCAA","TTCCTCAAGG","TACTGGATAA","GGATGGTCCG","TGCGCGGCAT","AGCGCGGAAC","CTGGATGCCA","TTATCTTATT","GTCTTGATGC","TGGTCAAGAT","CCGTCAGTTG","TGACCTTAAT","CTGATTCTCT","CCTTCCTTCA","CATACGGATT","GAGCTTAGAA","CTTATCGCAG","CCGGCAATCC","GTTCTATTGG","CAATCGACGG","GGAGGAGCTT","GGTTATATTA","GATGGTTACC","AGCGTCTGGT","ATGGTACTCT","CGTAGGTTGG","ATCGGTTCCG","ATAGGTTATC","CGAGGCTACC","GCCATAGTAG","GGCCGTTGAT","TTCGACCAGA","CGCCTAGAGA","GTTCTCCGGA","CGCCTCGGCG","TGGATCTGCA","GTCTATACCA","TTATCCATCC","TGGACCGCCA","ACTCGCCGCT","TAGCATTGAT","GTCATCGACC","GTATAACTCC","CGGCGAGTCC","AGCCTGACTG","ACGAGCCTCT","ATGAAGGAGG","TGCTTCCGCC","ATGCGCGAGT","TATGAGTTGA","ATGGCTCGGT","GCTGGCTGAA","CTGGATTGGA","GTTGGCATGG","AGAATGCTGG","GGAGAGATCC","GTCAACTAAT","GGTACCAACC","ACGGCAAGCA","CCGGATTACT","GAAGGTTCGG","AGTAAGTTAC","GGCCTAGCGG","ACTGATGAGT","GCGGTCAGAG","TAGGTATGAA","GAGGCAAGGC","TTCGAACGTC","GTCCTGCCGA","TAGGCGATAT","GAACCTCTAA","GGAACGCCTA","GGTATTAGCG","TGCAAGCGCG","CCAACTGAAT","GGCATATATA","AGCCGAGGAC","AATGCTCTGA","GTATTGAAGT","TTAATGACTT","GTATCTCGAT","AGCCGTTCCA","GTACGATACT","ATGGTTAGAA","CCTTGGCGAA","TTGGCTACTA","CCGGAGGACC","AGAGAGTTGG","TTACGCCATT","ATATGCCGCG","CCTGGATAGT","GAGGCGTAAT","CGACGGAGGC","ATAAGTCGAA","TCGAGGTCTT","TCAGACTGGA","GTCGATTATG","CAGGACCAGC","GCGCAGGCGA","GTTAGCTAAG","CTACCATCTC","CTTACTTCCT","GACTTCTACT","TAACCGCTGA","CTTGCCAATA","ATAGTATTGC","CCGCGGTCAA","ATCTCCTGAA","CGGACCGGTA","ATAACCGTCC","TTCTGATTAA","GCCTACGTTC","CATTGATACC","AAGTCTTCTG","AGTTATTGAC","ATTAGTAACG","GTCATATTCG","CGACTTCGGA","TATCGTAGGA","ACCGCTCAAC","AGCTGCGTTG","TATTATTGGT","GCTATGCGAT","TAGGCTTCGA","CGTATCTGCT","TGGTTGAGGA","GTTCCAACTT","CATGGCGTCA","GTCGTTGCTT","ATTATACTAC","CGAATAGTTG","GGTAATGGAC","GGAGTCCATT","GATGAACGGT","GCGGTCTTAT","CTCCGAATGC","TGAGAACTTG","ACTCCGACCA","TAATCCTCTT","CAGTATGGTC","ACGTCTCCTC","CTCCAGGCGG","AAGATTCCTC","TCTCGAAGGC","ACTAGTCTCG","TCCAGCTCTC","TAACCTGCAA","GTTCTGACGT","CCAAGGTCGG","CGTTACTTAG","CGGCATGGAC","GTCTGCGCTT","TGAATTATGA","AAGCCGGTAA","CGCATTAATT","GGCCGACCGT","AGCGAAGTAG","GGCCATACTT","GAGAGTTAAC","GTTGGCGACC","AGGTTATACC","TTCTATAGTT","GAGGTCTGAC","GAGAATGAAG","GTCCGCAGTA","TTGGTTAACC","TTATGGCCGG","CAGACTTCAT","GTCGCTGGAG","TCGAGAAGGT","TCGGCGTACC","CAACGCCTGG","GCTGCCATTC","AAGTTCGTAA","TCATGCTCAG","TATCTCTCAT","GAGGCTTGGT","ACCAATATGG","GGAACTGAGC","GCCGAACTGC","TAATACGGAC","GGATATACGG","TTACTCTCTT","CGCAACGCCA","GTATAAGGCA","AGATAATTCC"]
pcr_i7_list_384 = ["GAATGCCTTG","CGGACTGGCC","CTGGCAGCGG","TCATCAGCCA","TAACCAGTTA","AGCTCATTCG","CAACGGTCTT","TGGTATCAGA","GGTAACTACG","CTTCAGGCCA","ATCCTTCAAC","TGATTGAATG","TGCTTAGCGA","CGAGTTCGCC","CGGCATCTTG","TCAACGAGGT","GGTCTAGGAA","AAGTCAACGG","TAACCAATCG","ATCGGCTTGA","CAATCTCAAT","CAGTTGCTTG","ATGCTAACCT","CGGATCTCCG","TGACGCGACC","TCTGACGAAC","GATCATGATA","TTCTCTACGA","TACGCAGGTT","TGCGAATCGG","CTAAGAGTTA","TCCATACCTG","ACCGGTCTGC","TAGAACCTCC","TGACCATCAG","GACTCGAGGA","GCATTAGGCG","GGACTTATGG","TTGGTTCTCA","TAATTCCGGT","CTACTGCAAG","CGTAGCCGTA","TTGACTTCGT","TCAGCGCAAT","ATATGGTACC","TGCCGGTACC","GGTACTTCCA","AATTGAGTTA","GCGGCGCAAG","AATACCTGCC","AGCATACGGC","GCCTCGGTAA","ACGGACCTGC","CATTCTGATG","GCTTACTATA","TCCAATGGTT","ACGGTAATTG","TCGGAGTATG","GGTTAGTTGG","CTGCATTACG","TGCGGCCTGG","CATCTATTCG","GCAGCTGATT","ATATCTTCCG","GTCGTACGCG","TTGATGCCTC","TTGGAACTGG","GGAGAATGGC","TTCCATTCTT","CCGAAGGTTA","CCGGCTGAGC","ACGAGGAGCC","AAGATCAGTT","ATCCGAACGG","CTGCGCTGGT","TGCCATGGAA","AGCTGAACGC","GAGACGTTCT","CTGACTACTT","CGCCGATATC","TGGATCGTTA","GTACCTGAAT","GAATTCTCCA","AGACCAGGTT","CTAGCTTCTT","GGCAACCTTC","CCGTTAGCAA","GGCGTCTGCC","CGCGCAACCG","AGTCTTAATT","TAAGGAACGG","TCGTTCATTA","ACTCCAGATT","CGATATATCT","GGCTGGCTCT","GGCGGATACC"]
pcr_i5_list_384 = ["GACCTCCTTG","ACCAGTTCAG","TTCTGGCGCA","GGCAATCAAG","TTGGCCAGGT","TCTCCTCCTG","CTTCAGTAGT","AAGGTTAAGT","TGAGCCGCGG","TGGCCGTCTT","AACCTTCTCT","GGCGCCATCG","GTCGACCATT","CATTACCAAG","CTAGTTAAGT","CATATTGCAT","GCCAAGGCAA","GATCCATAAC","ATGAATACCA","TTATTGCTAA","TATTCCTTGC","GGTACCGGCT","GTACTGGTTG","CGTTGATGAT","ATGCCTCTTC","TTATATTGAA","GTTCGTCTGA","ATCCGGATAA","TTGCCTTATC","ATCTAGTTAG","AGATAGAGGT","TTCTAATCCG","AAGATGGATT","TTCCAAGTAA","ATGAAGCTTG","ATCCTTGGAG","GTTCCGTATC","GCCGGAGCGG","AGATGGCGCA","TTATTAAGAA","CTCTTCCAGG","TTAGCGGTTG","AGTAGATCTA","GTATCAATAT","GGTCTTCGGT","CCTAGTCTAT","GACCGTTACT","TATGCTTACT","CAGGAATCGG","AGGAGTAAGG","TGAGGTTAGC","GCATCGCGCA","GAACGTCAGG","TGAGTTGGAC","ATAGCATCAA","TCTAGTTCCA","TGCGTCTATA","ATAGAGAGCC","TCTTAATCAG","CGGAAGACTT","GCCAGACTCA","TAACTGGCTT","AAGGCCATCA","GTATACGTAA","CCAATTCCAT","AATATAGCGG","GTTCGAGGCG","ACGTACTAGG","ATCGAATCTG","TCTGCGTCGC","TAAGCTCGCC","GGCATAACCG","AGTCGCGTCG","TGGACCAAGG","CTGATCAGAG","TCTCATTGCC","ACTCGGCATA","GAGCTCAGCC","ATTGAGGATA","TTAGCGGAAC","CCTTCATCCG","GGCAGCAGTT","TAAGTAAGTA","TATCAAGGTC","GCGGTTGGAA","GAATGAATAA","GAAGAGTATT","TACCTTAGCT","CGTCAGTCAA","GCGCCTAGTT","GTAGTAGTCC","ACGCTTCGTC","GATAGCCGAT","AGTCATGGAT","TCTTAAGGAC","CCGGTCCTAA"]
lig_i5_list_384 = ["AAGTAGACTA","ATGGAAGCAT","TGGATCAGGC","GTTACTTAGC","ACCGCCGCAA","CTCAAGTCCT","CGGTCGACTA","TTCGCCGTAA","GCTCCGCTTG","ACTTAAGATA","GGCATGGCCA","CTTCGGTATA","GAGATTCGCC","CTAGGCCGTT","GGCCAACGAT","ACGGAACCTG","TGATTCTCGT","TTGCGTCAAC","GAATGCAACC","TGCGGTTCAG","TTGGCCAACC","TTGGTTAAGC","CTTAAGTTCG","GTCCTCAGAA","GGAGTCGTCT","GGTACCTCTA","GATCGCTGAG","AGAGTACTCC","TCATTCTATT","GTTACTACCA","CCAGCTCGCC","CGCCGGTATG","TTAATTCGTA","GAAGGCTCCA","GAGACGTACG","GAAGAGCCTC","CCGATGCATA","GTAATGGTAT","TTCTATCTCA","GCAGCAGCTA","TTGCTCGATT","CCTCATCGGC","ACTTCAGCAA","AGGTCATCCT","AACGCGTCAG","CTATGCTTAC","GTTGCCGTTC","GCTTACCGCC","TGGCAAGTCA","CATCGAAGGA","AGAATCCTCG","GCAATCGGTT","CCTAAGATTC","CCTGCGCGCG","ATCAGCGCGA","GTACGATTCT","TTACCTTGCA","CCGGCTCAGC","TTCTGCAAGA","ATATACGCTT","CTCAGCAACC","CAATTCTAGG","ATCAGTCTCG","AATCCGCAAC","CGGTTACCTT","ACGTTAAGAC","CTATCCAACC","ATAAGCGAAT","CTTATATCGG","ATATGACGAC","TTACCGCATA","ATTCATCGCC","AGAAGCAGAA","GTTCGTCGTT","CATGCTTCCA","TCGGTACCAG","TTGAGCCAAT","AGATGACTGA","ACGCTAGAAG","GTTCAATTGC","GGACCGTCAA","CATTAACGGA","TAAGCAGTCC","CCGGTCAGTT","ATAACGGACT","ACGAGAAGAT","ATCCTCTTAA","AATCCAATAA","CTAGCAGGAT","TGGTCTCGGA","CCGAGTACTA","GATGACGAAG","GGCAGTCTTC","AATACGAATA","ACCTAGGAGA","GAAGCGCCAA","CGTTACGTTG","GTCGCGAATA","TTAGAGCCTG","ACGGTCATCA","ACGTAGCAGG","CGACCGAGAG","AAGCGGTTCT","TCGGAATAAC","AAGTTCGCTG","AATAATCGGT","AGGCGAAGGC","AAGCCGCCGC","TCGGCCGATG","AGCGACTGCT","CTTAATGAGC","AATTCCTCTC","GCTGGTCTCC","AGTATTGCTA","TCTAGGATAA","GGTCCTGCAA","CGCTTCAATT","GGATTATTAT","TCCGGCTGAT","CCGCCTCGTT","TTATAATCAA","CCATTGAACG","ACTCCAACGG","ACCTCCTGAA","AGAGGCCGGC","CTGCCTCTTC","CAGTATCCTT","GTCAACTAGC","TGACGCAGTC","GTCAATACGA","TGAACTTCGA","CGTACCAACG","AGAGATGAAT","TATTCCAATT","GGATGCGATT","GTAACCAGGT","CCTCGTCATA","AATGGTCTTA","ATGAATGCCT","GTCCGTAGAT","CCATCCTAGT","TGGTTCCTAC","GCGCCTTCCG","CGTACTACGC","GGCCGCGGTT","TGGATAGTTG","CGGCGCCAGG","CAAGCTCAGG","ATCATCCTTC","CCTCCGGAGT","CCATTGCTGG","TATTCGCAGT","CCGGTTAAGT","ATATTCTACC","TACGGATCGT","TTCTCTCCAG","CCAAGAGCAA","TTGGTTCGAG","AACGGATTAC","CATCTTCAGA","TTGAACCTCC","GGAATTCCAA","ATAGGTCCAA","GCCATGGTAC","GAGCTCTTCA","CCGAGGCAAC","GTCTCTAGTT","GCTGGTTATA","TCGTAGGTCA","AACTCAGACG","TGCTGCCGGA","TGGAGGCAAG","ACTGATGCGA","ACGACTCCTC","TGGCAGCGAA","CCGATACTCT","CAATATAGGC","ACCGGCCGAC","AATAAGGCTC","CATCATAGCA","GATGATCCAT","ATGGCAATAC","ACCAGAACCA","GGTTCGACCT","CTTGGACGGA","CGGTCTCATA","AATCAGAGCC","CCTGAATACT","CTTGGAGACT","AAGACCTTAC","GCGAGCGCTC","TCGCAAGACG","CAATCTCGGA","TCGACCTACC","TTATAGGCAT","TTCTGGCCTA","AATTGGCGAT","AGGCAACCGC","ACGCAATGAT","ACTCGTTACT","CAAGTCAGCA","GGTACGGATA","CTCAATGGTC","ACTCCGAAGA","AATCCAAGGC","CGTCCATCGC","TTAGGTACGA","CTTACCATCT","CGAGATAAGA","CGATCTTCTA","GTTCCGATCC","GCCATAAGAA","AGATTCATAA","CAGAGCGTCA","GAGGAGCCAG","ATGAGGATCC","CGGCGCTCAA","CCTTACGTCA","TTCTCAGTCA","ATAGGAGCTT","TAGCCGACTC","AAGAACTCCG","CTTGAGACCG","AACGCCATAC","CCTACTCAAC","TCTAGCCTCC","CGCTTAGATC","GCTTAATCAT","AGCAGGTCAA","CCGCTTATAA","CCTAATGGAA","CATTGGTTCC","ACCAGTTATT","GACTCGCTTA","AATTGCTCCG","TGCGAATGCC","TGCCTCCAGC","TTAGGCTTGA","ACTCTGAGCT","TTCGAACGAC","TACTGCTCGG","TCTAGAAGAC","TATTGGAATA","AAGATATCAA","AGTCCATGGT","GGTTGAATAC","CTAGCGATGG","CCAATACGCC","GAGATACCTC","GCTCAGGAGC","ACTAGTTGCA","CCGCTATCCA","ATGCAACTTA","AGGACCAAGT","GGCGTCCTCA","GGAAGCTCGC","AATCGTTATA","CCGAGAGAGG","GAGGAACTCA","GGTCTGAGCA","ACCATTAAGC","AACTCTACCA","TGCTCAACTC","CCGGCAGCCT","CTCGCGACCA","TGGACCTTAG","GGAATGACCG","CAGTATGATG","CTCATGAATA","CTGCGTAGTA","GGAGCAACCT","TCAGTAATAC","GACGCTGCAT","ATCTCCAGCT","TTAGAACTGC","CCGTAACGCC","GACTATACCT","TGCGAGAGGT","ATAGGCCTCA","TCAATTCAAC","GGCTACCTGG","CTTCTTGAAC","TTACGCAGCA","ATAACCGCCA","CCAGTCGAAT","TAACCAACTA","TAGCTGGCGA","CAATGACTTG","CTCTGCGATC","GCATACCAAT","ACCTGATAGG","ATGGTTACCG","CTTAGGTTAT","CGCCTTAGCC","ACGTAACGTA","TTATTCTCTA","GCATAGTATA","CATGCGCATA","GTACCAATTC","AACGAGATCA","TCTCGATTAA","CGTCGACCGG","GTTGATAGCT","CAGCGTTGGT","ATCGGCATTG","GGTAGTCCTA","TAGCATCGCG","ACTGGTAACC","TCTACTGACC","CAAGAGTTAT","CAATTAGGAA","TTGCGGAACC","AGCCGAAGTA","GGTCCGTACG","AGCATGCAAG","TCTTCGTTGC","GTTGGAAGAA","GGAATACGGA","ACCGACGGTA","ATACTTGATT","GTCGTTCGCC","AATGATTGCA","CCATTAGCAG","CTGACTAATG","CCTTACGGAC","CGCATAACTA","AACCGGAGGA","AATGCAGCGG","CAGTCGGCAA","ATAGAATCCA","TCTCAGATAT","ACGAAGTATT","AGACTTATCA","TCGCGCCGTA","AGCTTGAAGA","CGGTAGCTAC","GCGGCATGCG","AAGACTGGCT","CGCATTCTTA","TACCGTCTCC","AATATAGTAT","ATCATAAGAT","ATCAATATCC","TATTACCAAC","GAGAAGACCA","TGGCGCTCTC","GCGACGATAA","AACGTCGCGA","ATGAAGCTTC","GCCATAGAGT","TGACTGAGAA","CGGCTCCGAA","CAGCCAATTG","GCGACCAGTT","CCTAACGACG","CTTCGCAATC","TGACTCCGTT","GATAGTCGCT","AAGGTACTAA","TTGCATGAGG","GTATTATATA","ATTCTTGGCT","TGCATCTTGG","GTTGGCTCAA","AATATCATTA","TAACTAAGTC","CAATAACCAA","TTATACTGCA","TGAGCAGAGC","TGCAAGCCAA","TGGAGAACGA","ATCGGATTCA","ACTAGACCGA","CGAGATGCTT","TTCTATTAAT","AATTAGTCCA","GCTCCAAGCC","CTCTTCCTAA","CCGCGTTAAC","GACGGAATAG","CTCAGAGTTG","GGACGTATGA","AAGATGAGTC","ATGCGCTACC"]

lig_i7_to_well_384 = get_well_dict(lig_i7_list_384, row_ordered=True)
pcr_to_well_384 = get_pcr_plate_dict(pcr_i5_list_384, pcr_i7_list_384)
lig_i5_to_well_384 = get_well_dict(lig_i5_list_384, row_ordered=False)

lig_i7_384 = set(lig_i7_list_384)
pcr_i7_384 = set(pcr_i7_list_384)
pcr_i5_384 = set(pcr_i5_list_384)
lig_i5_384 = set(lig_i5_list_384)

## Two-level Indexed TN5 barcode sets and well mappings
nex_i7_two_level_indexed_tn5_list = ["ATTACTCG", "TCCGGAGA", "CGCTCATT", "GAGATTCC", "ATTCAGAA", "GAATTCGT", "CTGAAGCT", "TAATGCGC", "CGGCTATG", "TCCGCGAA", "TCTCGCGC", "AGCGATAG"]
pcr_i7_two_level_indexed_tn5_list = ["TCGGATTCGG" ,"GCGGCTGCGG", "AGATTACGTT", "CTAACTAGGT", "CATAGCGACC", "CCGCTAAGAG", "ATGGAACGAA", "GCGTTCCGTT", "GGTTATCGAA", "GCATCGTATG", "AATACGATAA", "TTCCGTCGAC", "TCCGGCTTAT", "ACCAGGCGCA", "AGAGGAGAAT", "GTACTCCTAT", "GCTAACGGAT", "AGTTGAATCA", "TGATTAGGTA", "TCGTAGCATC", "TCTTGAGGTT", "AGGTCAGCTT", "TATTAGACTT", "CTCAATTAGT", "TCGCCGCCGG", "CCGTATGATT", "AACGCGCAGA", "CTCGTCGTAG", "CTAATTGCGA", "CGCGGCCATA", "AATATTACTT", "ATTGGCAGAT", "ATGGCGCCTG", "ATAAGGACTC", "TAGTAAGCCG", "ATTATGCAAG", "TTGGCAAGCC", "TTGATTGGCG", "GCATATGAGC", "GAACTCGACT", "CTAGCCAGCC", "TGCGACCTCT", "ATTCTTAGCT", "TTGATACGAT", "TATAATAGTT", "TTGCCGTAGG", "AGACCATATC", "TTGGTAAGGA", "CAGCTAGCGG", "CTAAGCCTTG", "CGTTACCGCT", "GACTGGACCA", "GCAAGACCGT", "TCAATCTCCT", "ATACCTCGAC", "TAGAGGCGTT", "TAGGTAACTT", "TTCGAATATT", "TGGACGACTA", "GTAGGCTGCA", "GTAGGATAAG", "CGTCGAGCGC", "ACTATTCATT", "TTGCTTAGAT", "CGAATGGAGC", "CTATATAGCC", "CTACTAATAA", "TGGTTGCCGT", "TCCTCTGCCG", "GATTCTTGAA", "GTAGCAGCTA", "CCTCAGCTCC", "AAGTAGCTCA", "TATTGCTGGA", "CCAGATACGG", "AACGAATTCG", "CGCTTATCGT", "AAGTACGCGA", "GATCTTCGCA", "TCTTAGCCTG", "TTATTGAGGC", "TTGCGAGCAT", "GCTTGAAGAG", "AGTCCGCTGC", "TAAGTCCTGA", "AGTTCTCATG", "CAGACTAAGG", "TCTATCGCTG", "GCGCTATGGT", "CATTATTATT", "AGCCGTAGTT", "TGATATTGCG", "ACGGCGTTAA", "GGCTTACTCC", "GCGCGTTCAT", "GAGCGCGATG"]
pcr_i5_two_level_indexed_tn5_list = ["CTCCATCGAG", "TTGGTAGTCG", "GGCCGTCAAC", "CCTAGACGAG", "TCGTTAGAGC", "CGTTCTATCA", "CGGAATCTAA", "ATGACTGATC", "TCAATATCGA", "GTAGACCTGG", "TTATGACCAA", "TTGGTCCGTT", "GGTACGTTAA", "CAATGAGTCC", "GATGCAGTTC", "CCATCGTTCC", "TTGAGAGAGT", "ACTGAGCGAC", "TGAGGAATCA", "CCTCCGACGG", "CATTGACGCT", "TCGTCCTTCG", "TGATACTCAA", "TTCTACCTCA", "TCGTCGGAAC", "ATCGAGATGA", "TAGACTAGTC", "GTCGAAGCAG", "AGGCGCTAGG", "AGATGCAACT", "AAGCCTACGA", "GTAGGCAATT", "GGAGGCGGCG", "CCAGTACTTG", "GGTCTCGCCG", "GGCGGAGGTC", "TAGTTCTAGA", "TTGGAGTTAG", "AGATCTTGGT", "GTAATGATCG", "CAGAGAGGTC", "TTAATTAGCC", "CTCTAACTCG", "TACGATCATC", "AGGCGAGAGC", "TCAAGATAGT", "TAATTGACCT", "CAGCCGGCTT", "AGAACCGGAG", "GAGATGCATG", "GATTACCGGA", "TCGTAACGGT", "TGGCGACGGA", "AGTCATAGCC", "GTCAAGTCCA", "ATTCGGAAGT", "GTCGGTAGTT", "AGGACGGACG", "CTCCTGGACC", "TAGCCTCGTT", "GGTTGAACGT", "AGGTCCTCGT", "GGAAGTTATA", "TGGTAATCCT", "AAGCTAGGTT", "TCCGCGGACT", "TGCGGATAGT", "TGGCAGCTCG", "TGCTACGGTC", "GCGCAATGAC", "CTTAATCTTG", "GGAGTTGCGT", "ACTCGTATCA", "GGTAATAATG", "TCCTTATAGA", "CCGACTCCAA", "GCCAAGCTTG", "CATATCCTAT", "ACCTACGCCA", "GGAATTCAGT", "TGGCGTAGAA", "ATTGCGGCCA", "TTCAGCTTGG", "CCATCTGGCA", "CTTATAAGTT", "GATTAGATGA", "TATAGGATCT", "AGCTTATAGG", "GTCTGCAATC", "CGCCTCTTAT", "GTTGGATCTT", "GCGATTGCAG", "TGCCAGTTGC", "CTTAGGTATC", "GAGACCTACC", "ATTGACCGAG"]
nex_i5_two_level_indexed_tn5_list = ["TATAGCCT", "ATAGAGGC", "CCTATCCT", "GGCTCTGA", "AGGCGAAG", "TAATCTTA", "CAGGACGT", "GTACTGAC"]

nex_two_level_indexed_tn5_list = get_row_col_matrix(nex_i5_two_level_indexed_tn5_list, nex_i7_two_level_indexed_tn5_list)
nex_two_level_indexed_tn5_to_well = get_well_dict(nex_two_level_indexed_tn5_list, row_ordered=True)
pcr_two_level_indexed_tn5_to_well = get_pcr_plate_dict(pcr_i5_two_level_indexed_tn5_list, pcr_i7_two_level_indexed_tn5_list)

nex_i7_two_level_indexed_tn5 = set(nex_i7_two_level_indexed_tn5_list)
pcr_i7_two_level_indexed_tn5 = set(pcr_i7_two_level_indexed_tn5_list)
pcr_i5_two_level_indexed_tn5 = set(pcr_i5_two_level_indexed_tn5_list)
nex_i5_two_level_indexed_tn5 = set(nex_i5_two_level_indexed_tn5_list)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A program to correct cell barcodes and split FASTQs by sample for sciATAC data.')
    parser.add_argument('-1', '--input1', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R1.')
    parser.add_argument('-2', '--input2', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R2.')
    parser.add_argument('--samplesheet', required=True, help='Samplesheet describing the layout of the samples.')
    parser.add_argument('--output1_prefix', default='r1.', help='Prefix to add to R1 output files. Must be different from R2 prefix.')
    parser.add_argument('--output2_prefix', default='r2.', help='Prefix to add to R2 output files. Must be different from R1 prefix.')
    parser.add_argument('--stats_out', required=True, help='JSON file with output stats about processed reads and correction rates.')
    parser.add_argument('--two_level_indexed_tn5', action='store_true', help='Flag to run assuming that the library is a two-level indexed-TN5 sample.')
    parser.add_argument('--wells_384', action='store_true', help='Flag to run assuming That the known barcode set is the 384 well set.')
    parser.add_argument('--well_ids', action='store_true', help='Flag to output cell IDs that are composed of well IDs rather than the actual sequences.')
    parser.add_argument('-X', '--nextseq', help='Specify this flag if the run was on NextSeq or similar (requires reverse complementing of some barcode sequences). Not required for NovaSeq, for example.', dest='nextseq', action="store_true")
    args = parser.parse_args()

    if args.output1_prefix == args.output2_prefix:
        raise ValueError('output1_prefix and output2_prefix arguments must not be the same.')

    if args.two_level_indexed_tn5 and args.wells_384:
        raise ValueError('There is no 384 well barcode set for indexed Tn5, may not specify both --two_level_indexed_tn5 and --wells_384.')

	# Set up the right index set depending on the indices
    if args.two_level_indexed_tn5:
        tagi7 = nex_i7_two_level_indexed_tn5_list
        pcri7 = pcr_i7_two_level_indexed_tn5_list
        pcri5 = pcr_i5_two_level_indexed_tn5_list
        tagi5 = nex_i5_two_level_indexed_tn5_list
    else:
        if args.wells_384:
            tagi7 = lig_i7_list_384
            pcri7 = pcr_i7_list_384
            pcri5 = pcr_i5_list_384
            tagi5 = lig_i5_list_384
            lig_i7_to_well = lig_i7_to_well_384
            lig_i5_to_well = lig_i5_to_well_384
            pcr_to_well = pcr_to_well_384
        else:
            tagi7 = lig_i7_list
            pcri7 = pcr_i7_list
            pcri5 = pcr_i5_list
            tagi5 = lig_i5_list
            lig_i7_to_well = lig_i7_to_well
            lig_i5_to_well = lig_i5_to_well
            pcr_to_well = pcr_to_well

    # Build up sample mapping from indices to samples
    index_mask, sample_lookup = get_sample_lookup(open(args.samplesheet), pcri7, tagi7, tagi5, pcri5)

    tagmentation_i7_whitelist = tagi7
    pcr_i7_whitelist = pcri7

    if not args.two_level_indexed_tn5:

        if args.nextseq:
            p5_pcr_rc_map = {reverse_complement(k):k for k in pcri5}
            p5_tagmentation_rc_map = {reverse_complement(k):k for k in tagi5}

            pcr_i5_whitelist = set([reverse_complement(x) for x in pcri5])
            tagmentation_i5_whitelist = set([reverse_complement(x) for x in tagi5])
        else:
            pcr_i5_whitelist = pcri5
            tagmentation_i5_whitelist = tagi5

    else:
        tagmentation_i7_whitelist = nex_i7_two_level_indexed_tn5
        pcr_i7_whitelist = pcr_i7_two_level_indexed_tn5

        if args.nextseq:
            p5_pcr_rc_map = {reverse_complement(k):k for k in pcr_i5_two_level_indexed_tn5}
            p5_tagmentation_rc_map = {reverse_complement(k):k for k in nex_i5_two_level_indexed_tn5}

            pcr_i5_whitelist = set([reverse_complement(x) for x in pcr_i5_two_level_indexed_tn5])
            tagmentation_i5_whitelist = set([reverse_complement(x) for x in nex_i5_two_level_indexed_tn5])
        else:
            pcr_i5_whitelist = pcr_i5_two_level_indexed_tn5
            tagmentation_i5_whitelist = nex_i5_two_level_indexed_tn5

    tagmentation_i7_correction_map = construct_mismatch_to_whitelist_map(tagmentation_i7_whitelist, 2)
    pcr_i7_correction_map = construct_mismatch_to_whitelist_map(pcr_i7_whitelist, 2)
    pcr_i5_correction_map = construct_mismatch_to_whitelist_map(pcr_i5_whitelist, 2)
    tagmentation_i5_correction_map = construct_mismatch_to_whitelist_map(tagmentation_i5_whitelist, 2)

    # Set up all input/output files
    output_files = {}
    for sample in list(set(sample_lookup.values())):
        output_file_1 = '%s%s.fastq' % (args.output1_prefix, sample)
        output_file_2 = '%s%s.fastq' % (args.output2_prefix, sample)
        output_files[sample] = {}
        output_files[sample]['r1'] = open(output_file_1, 'w')
        output_files[sample]['r1_name'] = output_file_1
        output_files[sample]['r2'] = open(output_file_2, 'w')
        output_files[sample]['r2_name'] = output_file_2

    if1 = FastqGeneralIterator(args.input1)
    if2 = FastqGeneralIterator(args.input2)

    totreads = 0
    total_not_specified_in_samplesheet = 0
    validreads = {}
    validreads['pcr_i5'] = 0
    validreads['pcr_i7'] = 0
    validreads['tagmentation_i5'] = 0
    validreads['tagmentation_i7'] = 0
    validreads['all_barcodes'] = 0

    start = time.time()

    for (r1_name, r1_seq, r1_qual),(r2_name, r2_seq, r2_qual) in zip(if1, if2):
        totreads += 1

        # Get barcodes and correct
        tagmentation_i7_seq, pcr_i7_seq, pcr_i5_seq, tagmentation_i5_seq = get_barcode_seqs(r1_name, args.nextseq, args.two_level_indexed_tn5)
        tagmentation_i7_seq = correct_barcode(tagmentation_i7_seq, tagmentation_i7_correction_map)
        pcr_i7_seq = correct_barcode(pcr_i7_seq, pcr_i7_correction_map)
        pcr_i5_seq = correct_barcode(pcr_i5_seq, pcr_i5_correction_map)
        tagmentation_i5_seq = correct_barcode(tagmentation_i5_seq, tagmentation_i5_correction_map)

        # Skip invalid reads and track valid read count for error checking
        if tagmentation_i7_seq is not None:
            validreads['tagmentation_i7'] += 1
        if tagmentation_i5_seq is not None:
            validreads['tagmentation_i5'] += 1
        if pcr_i7_seq is not None:
            validreads['pcr_i7'] += 1
        if pcr_i5_seq is not None:
            validreads['pcr_i5'] += 1
        
        if tagmentation_i7_seq is None or pcr_i7_seq is None or pcr_i5_seq is None or tagmentation_i5_seq is None:
            continue

        validreads['all_barcodes'] += 1

        # Map back to original whitelist if on nextseq so barcodes are always same on every sequencer
        if args.nextseq:
            pcr_i5_seq = p5_pcr_rc_map[pcr_i5_seq]
            tagmentation_i5_seq = p5_tagmentation_rc_map[tagmentation_i5_seq]

        # Convert to well IDs if requested
        # Note that for two level Tn5 barcode well comes first then PCR,
        # for three-level will be Tn5 N7, Tn5 N5, PCR WELL ID
        if args.two_level_indexed_tn5:
            barcodes_string = get_two_level_barcode_string(tagmentation_i7_seq, pcr_i7_seq, pcr_i5_seq, tagmentation_i5_seq, nex_two_level_indexed_tn5_to_well, pcr_two_level_indexed_tn5_to_well, args.well_ids)
        else:
            barcodes_string = get_barcode_string(tagmentation_i7_seq, pcr_i7_seq, pcr_i5_seq, tagmentation_i5_seq, lig_i7_to_well, lig_i5_to_well, pcr_to_well, args.well_ids)

        sample_index = tuple([index for use_index,index in zip(index_mask, [pcr_i7_seq, tagmentation_i7_seq, tagmentation_i5_seq, pcr_i5_seq]) if use_index])
        sample = sample_lookup.get(sample_index, None)

        if not sample:
            total_not_specified_in_samplesheet += 1
        else:
            output_files[sample]['r1'].write(''.join(['@', barcodes_string, ':', str(totreads), '\n', r1_seq, '\n+\n', r1_qual, '\n']))
            output_files[sample]['r2'].write(''.join(['@', barcodes_string, ':', str(totreads), '\n', r2_seq, '\n+\n', r2_qual, '\n']))

    if totreads == 0:
        raise ValueError('No reads found in fastq input.')

    # Output basic stats
    for stat in validreads:
        validreads[stat] = validreads[stat] / float(totreads)
    validreads['total_input_reads'] = totreads
    validreads['total_not_specified_in_samplesheet'] = total_not_specified_in_samplesheet

    with open(args.stats_out, 'w') as f:
        f.write(json.dumps(validreads, indent=4))

    # Error checking and compress output
    if validreads['all_barcodes'] < 0.05:
        raise ValueError('Warning, you had less than 5 percent of all reads pass index correction. Something may have gone wrong here w.r.t. index sets or the expected library configuration not matching the data...')

    print('Done correcting barcodes in %s minutes. Starting compression...' % ((time.time() - start) / 60.0))
    start = time.time()
    for sample in output_files:
        output_files[sample]['r1'].close()
        output_files[sample]['r2'].close()
        subprocess.check_call('gzip -f %s' % output_files[sample]['r1_name'], shell=True)
        subprocess.check_call('gzip -f %s' % output_files[sample]['r2_name'], shell=True)
    print('Done compressing with gzip in %s minutes.' % ((time.time() - start) / 60.0))
