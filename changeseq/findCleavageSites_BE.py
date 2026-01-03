import argparse
import HTSeq
import os
import pyfaidx
import regex
from statsmodels.distributions.empirical_distribution import ECDF
import numpy as np
import pysam
import pandas as pd
import traceback
import logging
import pickle

logger = logging.getLogger('root')
logger.propagate = False


# reference_start is 0-index, reference_end is 1-index
def is_outie(x, y):
    ## simply the overlap
    if x.is_reverse:
        if not y.is_reverse:
            if x.reference_start < y.reference_start:
                return True, x.reference_end - y.reference_start
    else:
        if y.is_reverse:
            if y.reference_start < x.reference_start:
                return True, y.reference_end - x.reference_start
    return False, 0


# make sure BWA version >= 0.7.7
def get_conversion_pos(r, overlap_length=7):
    if overlap_length < 0:
        return []
    out = []
    actual_start = 0
    # expected pos is len(overlap) + 2, 1-index, we gave it 2bp flank
    expected_min = overlap_length
    expected_max = overlap_length + 4
    try:
        myList = r.get_aligned_pairs(
            with_seq=True)  # [(0, 82316526, 'G'), (1, 82316527, 'G'), (2, 82316528, 'C'), (3, 82316529, 'T'), (4, 82316530, 'G'),
    except Exception as e:
        logger.error(traceback.format_exc())
        return []
    # qualities = r.query_qualities
    if r.is_reverse:
        myList = myList[::-1]
    q_seq = r.query_sequence  # qseq == 'GGCTGCTGGGGGCAGTGGGCTGGGAGCT
    for q, ref_pos, base in myList:
        if q == None:
            continue
        actual_start += 1
        # if qualities[q]<=10:
        # 	continue
        if ref_pos == None:  ## this is chrom position
            continue
        if expected_min <= actual_start <= expected_max:
            if q_seq[q] == "G" and base.upper() == "A":  ### does the query seq, so reference, match the base
                out.append(ref_pos)
            if q_seq[q] == "C" and base.upper() == "T":
                out.append(ref_pos)
        if actual_start > expected_max:
            break
    return out


def get_read_strand(read):
    if read.is_reverse:
        return "-"
    return "+"


def get_read_start(read):
    if read.is_reverse:
        return read.reference_end
    return read.reference_start + 1


def tabulate_start_positions_BE(bam=None, label=None, output_dir=None,  # inputs
                                BEmodel_min_overlap=2, BEmodel_max_overlap=15, mapq_threshold=0):  # other filters

    # print (bam)
    # output files
    SE_read_stats_file = f"{output_dir}/{label}.SE.read_stats.csv"
    SE_read_stats_list_of_dict = []  # # read, chr, pos, converted_pos_list
    PE_read_stats_file = f"{output_dir}/{label}.PE.read_stats.csv"
    PE_read_stats_list_of_dict = []  # # read, chr, pos1, pos2, is_outie, overlap_bp, converted_pos_list

    # expected conversion position is overlap_bp+2

    # housekeeping
    logger.info("Reading bam to a dictionary")
    reads_dict = bam_to_dict(pysam.AlignmentFile(bam, "rb",threads=48), mapq_threshold)
    n_reads = len(reads_dict.keys())
    logger.info("Finished bam to dictionary")
    ga_coverage_paired = HTSeq.GenomicArray("auto", stranded=False, typecode="d")  # same as old ga_coverage
    ga_coverage_single = HTSeq.GenomicArray("auto", stranded=False,
                                            typecode="d")  # subset of ga_noise, save this var in case we want to rescue some reads
    ga_converted = HTSeq.GenomicArray("auto", stranded=False, typecode="d")
    ga_noise = HTSeq.GenomicArray("auto", stranded=False, typecode="d")

    # main
    read_count = 0
    not_noise_count = 0
    PE_read_count = 0
    for read in reads_dict:
        read_count += 1
        if not read_count % 150000:
            logger.info(f'processing {read_count} ; {round(read_count/ n_reads,4)*100}%')
        read1_list = reads_dict[read][0]
        read2_list = reads_dict[read][1]

        # parse single-end reads
        if len(read1_list) == 0 or len(read2_list) == 0:
            # if len(read1_list) == 0:
            #     ending = "/1"
            # if len(read2_list) == 0:
            #     ending = "/2"
            # SE_read_count += 1
            # read_list = read1_list + read2_list
            # for r in read_list:
            #     read_start = get_read_start(r)
            #     # ga_coverage_single[HTSeq.GenomicPosition(r.reference_name, read_start)]+=1
            #     # converted_pos_list = get_conversion_pos(r)
            #     converted_pos_list = []
            #     SE_read_stats_list_of_dict.append({"read": read + ending,
            #                                        "chr": r.reference_name,
            #                                        "pos": read_start,
            #                                        "strand": get_read_strand(r),
            #                                        "converted_pos_list": converted_pos_list})
            # # print (f"{read} (SE) is noise.")
            # # ga_noise[HTSeq.GenomicPosition(r.reference_name, read_start)]+=1
            # # ga_noise[HTSeq.GenomicInterval(r.reference_name, r.reference_start, r.reference_end)] += 1

            pass
        else:
            PE_read_count += 1
            noise_flag = True  # only count noise once
            # enumerate all possible pairs
            for i in read1_list:
                read1_start = get_read_start(i)
                read1_strand = get_read_strand(i)
                for j in read2_list:
                    read2_start = get_read_start(j)
                    read2_strand = get_read_strand(j)

                    # nonmatching chromosomesq
                    if i.reference_name == j.reference_name:

                        flag, overlap_bp = is_outie(i, j)
                        # outie and within overlap distance
                        if flag:
                            R1_converted_list = get_conversion_pos(i, overlap_bp)
                            R2_converted_list = get_conversion_pos(j, overlap_bp)
                            converted_pos_list = R1_converted_list + R2_converted_list

                            if BEmodel_min_overlap<= overlap_bp <= BEmodel_max_overlap:
                                noise_flag = False
                                ga_coverage_paired[HTSeq.GenomicPosition(i.reference_name, read1_start)] += 1
                                ga_coverage_paired[HTSeq.GenomicPosition(i.reference_name, read2_start)] += 1
                                for pos in converted_pos_list:
                                    ga_converted[HTSeq.GenomicPosition(i.reference_name, pos)] += 1
                            else:
                                noise_flag = True
                        else:
                            converted_pos_list = []
                            noise_flag = True

            if noise_flag:
                ga_noise[HTSeq.GenomicPosition(i.reference_name, read1_start)]+=1
                ga_noise[HTSeq.GenomicPosition(i.reference_name, read2_start)]+=1

            else: #moved this to be only stoed when overlap and outie are present
                # read, chr, pos1, pos2, is_outie, overlap_bp, converted_pos_list
                not_noise_count +=1
                PE_read_stats_list_of_dict.append({"read": read,
                                                   "chr": i.reference_name,
                                                   "pos1": i.reference_start,
                                                   "strand1": read1_strand,
                                                   "pos2": j.reference_start,
                                                   "strand2": read2_strand,
                                                   "is_outie": flag,
                                                   "overlap_bp": overlap_bp,
                                                   "converted_pos_list": converted_pos_list})

    pd.DataFrame.from_dict(PE_read_stats_list_of_dict).to_csv(PE_read_stats_file, index=False)

    logger.info(f"{len(reads_dict.keys())} reads were mapped and met MAPQ threshold: {mapq_threshold}")
    logger.info(
        f"{len(PE_read_stats_list_of_dict)}, {round(len(PE_read_stats_list_of_dict) / n_reads, 4) * 100}% outie paired reads were > {BEmodel_min_overlap} and < {BEmodel_max_overlap}")
    print(f"{not_noise_count}, {round( not_noise_count/ n_reads, 4) * 100}% outie paired reads were > {BEmodel_min_overlap} and < {BEmodel_max_overlap}")

    return ga_coverage_paired, ga_coverage_single, ga_converted, ga_noise, PE_read_count,not_noise_count


""" Find genomic windows (coordinate positions)
"""


def find_windows(ga_windows, window_size):
    # Initialize comparison position
    last = HTSeq.GenomicInterval("0", 0, 0)
    # Iterate through window GenomicArray and consolidate windows that are within 3 bp, up to a maximum of 10 bp.
    for iv, value in ga_windows.steps():
        if value:
            if iv.chrom != last.chrom or iv.start - last.end > window_size or iv.end - last.start > 10:
                last = iv
            else:
                consolidated_interval = HTSeq.GenomicInterval(iv.chrom, last.start, iv.end)
                ga_windows[consolidated_interval] = 1
                last = consolidated_interval

    return ga_windows  # Return consolidated GenomicArray


""" Find actual sequences of potential off-target sites
"""


def output_alignments(narrow_ga,
                      ga_windows,
                      narrow_ga_converted,
                      ga_narrow_windows_noise,
                      control_narrow_ga,
                      control_narrow_ga_converted,
                      control_ga_narrow_windows_noise,
                      reference_genome,
                      target_sequence,
                      target_name,
                      bam_filename,
                      edist_threshold,
                      mismatch_threshold,
                      bulge_threshold,
                      BEsearch_radius,
                      pkl_data,
                      out):


    pkl_out = '{0}_total_counts.pkl'.format(out)
    outfile_matched = '{0}_identified_matched.txt'.format(out)
    outfile_unmatched = '{0}_identified_unmatched.txt'.format(out)
    # dictionary to store the matched reads
    matched_dict = {}
    # dictionary to add read_count for each pair chromosome:start_position among matched reads
    reads_dict = {}
    # dictionary to store window_start. For duplicated matched off-target.
    control_reads_dict = {}
    # dictionary to store window_start. For duplicated matched control off-target.
    window_min = {}
    # dictionary to store window_end. For duplicated matched off-target.
    window_max = {}
    # dictionary to store the unmatched reads
    unmatched_dict = {}
    # radius is in additio to windows
    for iv, value in ga_windows.steps():
        if value:
            window_sequence = get_sequence(reference_genome, iv.chrom, iv.start - BEsearch_radius, iv.end + BEsearch_radius)

            offtarget_sequence_no_bulge, mismatches, offtarget_sequence_length, chosen_alignment_strand_m, start_no_bulge, end_no_bulge, \
                realigned_target, \
                bulged_offtarget_sequence, length, score, substitutions, insertions, deletions, chosen_alignment_strand_b, bulged_start, bulged_end = \
                alignSequences(target_sequence, window_sequence,edist=edist_threshold, bulges=bulge_threshold,max_score=mismatch_threshold)

            # get genomic coordinates of sequences
            mm_start, mm_end, b_start, b_end = '', '', '', ''
            if offtarget_sequence_no_bulge and chosen_alignment_strand_m == '+':
                mm_start = iv.start - BEsearch_radius + int(start_no_bulge)
                mm_end = iv.start - BEsearch_radius + int(end_no_bulge)
            if offtarget_sequence_no_bulge and chosen_alignment_strand_m == '-':
                mm_start = iv.end + BEsearch_radius - int(end_no_bulge)
                mm_end = iv.end + BEsearch_radius - int(start_no_bulge)

            if bulged_offtarget_sequence and chosen_alignment_strand_b == '+':
                b_start = iv.start - BEsearch_radius + int(bulged_start)
                b_end = iv.start - BEsearch_radius + int(bulged_end)
            if bulged_offtarget_sequence and chosen_alignment_strand_b == '-':
                b_start = iv.end + BEsearch_radius - int(bulged_end)
                b_end = iv.end + BEsearch_radius - int(bulged_start)

            #  define overall start, end and strand. For bed annotation purposes
            if offtarget_sequence_no_bulge:
                target_start_absolute = mm_start
                target_end_absolute = mm_end
                target_strand_absolute = chosen_alignment_strand_m
            elif not offtarget_sequence_no_bulge and bulged_offtarget_sequence:
                target_start_absolute = b_start
                target_end_absolute = b_end
                target_strand_absolute = chosen_alignment_strand_b
            else:
                target_start_absolute = iv.start
                target_end_absolute = iv.end
                target_strand_absolute = '*'

            name = iv.chrom + ':' + str(target_start_absolute) + '-' + str(target_end_absolute)
            read_count = int(max(set(narrow_ga[iv])))
            control_read_count = int(max(set(control_narrow_ga[iv])))
            filename = os.path.basename(bam_filename)
            full_name = str(target_name) +'_' + str(name)
            nuclease_converted_count =  str(int(max(narrow_ga_converted[iv])))
            control_converted_count =str(int(max(control_narrow_ga_converted[iv])))
            nuclease_noise_count = str(int(max(ga_narrow_windows_noise[iv])))
            control_noise_count = str(int(max(control_ga_narrow_windows_noise[iv])))


            if offtarget_sequence_no_bulge or bulged_offtarget_sequence:
                tag = iv.chrom + ':' + str(target_start_absolute)

                if tag not in reads_dict.keys():
                    reads_dict[tag] = read_count
                    control_reads_dict[tag] = control_read_count
                    window_min[tag] = [iv.start]
                    window_max[tag] = [iv.end]
                    matched_dict[tag] = [iv.chrom, target_start_absolute, target_end_absolute, name,target_strand_absolute,
                                         read_count,control_read_count,
                                         iv.start, iv.end, iv, window_sequence,
                                         offtarget_sequence_no_bulge, mismatches,
                                         chosen_alignment_strand_m, mm_start, mm_end,
                                         bulged_offtarget_sequence, length, score, substitutions, insertions, deletions,
                                         chosen_alignment_strand_b, b_start, b_end,
                                         filename, target_name, target_name, full_name, target_sequence,
                                         realigned_target,nuclease_converted_count, control_converted_count,nuclease_noise_count,control_noise_count]
                else:
                    current_read_count = reads_dict[tag]
                    current_control_read_count = control_reads_dict[tag]
                    reads_dict[tag] = max(current_read_count, read_count)
                    control_reads_dict[tag] = max(current_control_read_count,control_read_count)
                    window_min[tag].append(iv.start)
                    window_max[tag].append(iv.end)
                    matched_dict[tag] = [iv.chrom, target_start_absolute, target_end_absolute, name,target_strand_absolute,
                                         reads_dict[tag],control_reads_dict[tag],
                                         min(window_min[tag]), max(window_max[tag]), iv, window_sequence,
                                         offtarget_sequence_no_bulge, mismatches,
                                         chosen_alignment_strand_m, mm_start, mm_end,
                                         bulged_offtarget_sequence, length, score, substitutions, insertions, deletions,
                                         chosen_alignment_strand_b, b_start, b_end,
                                         filename, target_name, target_name, full_name, target_sequence,
                                         realigned_target, nuclease_converted_count, control_converted_count,nuclease_noise_count,control_noise_count]

            # no matching site but meets other critera
            else:
                untag = iv.chrom + ':' + str(iv.start)
                unmatched_dict[untag] = [iv.chrom, target_start_absolute, target_end_absolute, name, target_strand_absolute,
                                         read_count,control_read_count,
                                         iv.start, iv.end, iv, window_sequence,
                                         offtarget_sequence_no_bulge, mismatches,
                                         chosen_alignment_strand_m, mm_start, mm_end,
                                         bulged_offtarget_sequence, length, score, substitutions, insertions, deletions,
                                         chosen_alignment_strand_b, b_start, b_end,
                                         filename, target_name, target_name, full_name, target_sequence, 'none',
                                         nuclease_converted_count, control_converted_count,nuclease_noise_count,control_noise_count]

    # Write matched table

    # Yichao, add control reads
    logger.info("Writing matched table")
    tags_sorted = matched_dict.keys()
    tags_sorted = sorted(tags_sorted)

    pkl_data['total nuclease matched'] = [x for x in matched_dict.values()]
    pkl_data['total nuclease matched'] = len(tags_sorted)

    o1 = open(outfile_matched, 'w')
    # print('Chromosome', 'Start', 'End', 'Name', 'ReadCount', 'Strand',  # 0:5
    # 'MappingPositionStart', 'MappingPositionEnd', 'WindowName', 'WindowSequence',  # 6:9
    # 'Site_SubstitutionsOnly.Sequence', 'Site_SubstitutionsOnly.NumSubstitutions',  # 10:11
    # 'Site_SubstitutionsOnly.Strand', 'Site_SubstitutionsOnly.Start', 'Site_SubstitutionsOnly.End',  # 12:14
    # 'Site_GapsAllowed.Sequence', 'Site_GapsAllowed.Length', 'Site_GapsAllowed.Score',  # 15:17
    # 'Site_GapsAllowed.Substitutions', 'Site_GapsAllowed.Insertions', 'Site_GapsAllowed.Deletions',  # 18:20
    # 'Site_GapsAllowed.Strand', 'Site_GapsAllowed.Start', 'Site_GapsAllowed.End',  #21:23
    # 'FileName', 'Cell', 'Targetsite', 'FullName', 'TargetSequence', 'RealignedTargetSequence',  # 24:29
    # 'Position.Pvalue', 'Narrow.Pvalue', 'Position.Control.Pvalue', 'Narrow.Control.Pvalue','control_position_counts','control_window_counts',  # 30:33
    # sep='\t', file=o1)
    # Yichao Redefine output
    print('#Chromosome', 'Start', 'End', 'Genomic Coordinate', 'Strand','Nuclease_Read_Count', # 0:5 bed6 format
          'Control_Read_Count', 'Site_Sequence', 'Site_Sequence_Gaps_Allowed','Site_Substitution_Number', 'RNA_Bulge',
          'DNA_Bulge',  # contron window count, # 10:11, 15
          'File_Name', 'Cell', 'Target_site', 'Full_Name', 'Target_Sequence', 'Realigned_Target_Sequence',  # 24:29
          'Nuclease_Base_Converted', 'Control_Base_Converted','Nuclease_Noise', 'Control_Noise',
          'MappingPositionStart', 'MappingPositionEnd',
          'WindowName', 'WindowSequence',  # which column, -2 -1
          sep='\t', file=o1)
    o1.close()

    with open(outfile_matched, 'a') as o1:
        for key in tags_sorted:
            row = matched_dict[key]
            # '#Chromosome', 'Start', 'End', 'Genomic Coordinate', 'Strand','Nuclease_Read_Count', Control_Read_Count
            outline = [row[row_index] for row_index in [0, 1, 2, 3, 4, 5, 6]]
            # 'Site_Sequence', 'Site_Sequence_Gaps_Allowed',
            # 'Site_Substitution_Number','RNA_Bulge','DNA_Bulge',
            # 'File_Name', 'Cell', 'Target_site', 'Full_Name',
            # 'Target_Sequence', 'Realigned_Target_Sequence',  # 24:29
            #
            outline += [row[row_index] for row_index in
                        [11, 16,
                         19, 20, 21,
                         25, 26, 27, 28,
                         29, 30,
                         -4, -3, -2, -1]]
            # 'MappingPositionStart', 'MappingPositionEnd','WindowName', 'WindowSequence',  # which column, -2 -1
            outline += [row[row_index] for row_index in [7, 8, 9, 10]]
            pkl_data['total nuclease matched'] += int(row[5])
            pkl_data['total control matched'] += int(row[6])
            print(*(outline), sep='\t', file=o1)

    # Write unmatched table
    # print("Writing unmatched table", file=sys.stderr)
    untags_sorted = unmatched_dict.keys()
    untags_sorted = sorted(untags_sorted)

    with open(outfile_unmatched, 'w') as o2:
        for unkey in untags_sorted:
            unrow = unmatched_dict[unkey]
            pkl_data['total nuclease unmatched'] += int(unrow[5])
            pkl_data['total control unmatched'] += int(unrow[6])
            print(*(unrow), sep='\t', file=o2)


    # Save it as a pickle file
    with open(pkl_out, 'wb') as f:
        pickle.dump(pkl_data, f)


""" Reverse complement DNA sequence
"""


def reverseComplement(seq):
    compl = dict({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'X': 'X',
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n', 'x': 'x',
                  'Y': 'R', 'R': 'Y', 'W': 'W', 'S': 'S', 'K': 'M', 'M': 'K',
                  'y': 'r', 'r': 'y', 'w': 'w', 's': 's', 'k': 'm', 'm': 'k',
                  'D': 'H', 'H': 'D', 'V': 'B', 'B': 'V',
                  'd': 'h', 'h': 'd', 'v': 'b', 'b': 'v',
                  '.': '.', '-': '-', '_': '_'})
    out_list = [compl[bp] for bp in seq]
    return ''.join(out_list[::-1])


def regexFromSequence(seq,indels=1, errors=7):
    seq = seq.upper()
    """
    Given a sequence with ambiguous base characters, returns a regex that matches for
    the explicit (unambiguous) base characters
    """
    IUPAC_notation_regex = {'N': '[ATCGN]',
                            'Y': '[CTY]',
                            'R': '[AGR]',
                            'W': '[ATW]',
                            'S': '[CGS]',
                            'H': '[ACTH]',
                            'M': '[ACM]',
                            'B': '[CGTB]',
                            'K': '[GTK]',
                            'D': '[AGTD]',
                            'V': '[ACGV]',
                            'A': 'A',
                            'T': 'T',
                            'C': 'C',
                            'G': 'G'}

    pattern = ''
    for c in seq:
        pattern += IUPAC_notation_regex[c]

    pattern = "(" + pattern + "){s<=" + str(errors)+ ",i<=" + str(indels) + ",d<="+str(indels) +"}"
    return pattern

"""
Allow for '-' in our search, but do not allow insertions or deletions. 
"""


def extendedPattern(seq, errors):
    IUPAC_notation_regex_extended = {'N': '[ATCGN]', '-': '[ATCGN]', 'Y': '[CTY]', 'R': '[AGR]', 'W': '[ATW]',
                                     'S': '[CGS]', 'A': 'A', 'T': 'T', 'C': 'C', 'G': 'G'}
    realign_pattern = ''
    for c in seq:
        realign_pattern += IUPAC_notation_regex_extended[c]
    return '(?b:' + realign_pattern + ')' + '{{s<={0}}}'.format(errors)


"""
Recreate A!!! sequence in the window_sequence that matches the conditions given for the fuzzy regex. 
Currently only working for off-targets with at most one bulge !!! 
"""


## Yichao: Fix the code to get the 'realigned target sequence' to use the new feature in regex library that provides fuzzy_changes
def realignedSequences(targetsite_sequence, chosen_alignment, errors):
    m = chosen_alignment.group()
    q = targetsite_sequence
    substitutions, insertions, deletions = chosen_alignment.fuzzy_changes

    start = chosen_alignment.span()[0]
    indels = {}
    # deletion index is for targetsite_sequence
    # insertion index is for match sequence
    insertions = [i - start for i in insertions]
    deletions = [i - start for i in deletions]
    q_list = []
    m_list = []
    for i in insertions:
        indels[i] = True
    for d in deletions:
        indels[d] = False

    count = 0
    q_index = 0
    m_index = 0
    count = 0
    while q_index < len(q) or m_index < len(m):
        if q_index in indels:
            if not indels[q_index]:
                m_list.append("-")
                q_list.append(q[q_index])
                q_index += 1
                count += 1
                continue
        if m_index in indels:
            if indels[m_index]:
                q_list.append("-")
                m_list.append(m[m_index])

                ## update indel dict, I found the deletion positions are relative if insertion occurs first
                update_indels = []
                # for k in indels.keys():
                for k in list(indels):
                    if k > m_index:
                        if not indels[k]:
                            del indels[k]
                            update_indels.append(k - 1)
                for u in update_indels:
                    indels[u] = False
                m_index += 1
                count += 1
                continue
        m_list.append(m[m_index])
        q_list.append(q[q_index])
        q_index += 1
        m_index += 1
        count += 1

    realigned_target_sequence = "".join(q_list)
    realigned_offtarget_sequence = "".join(m_list)

    return realigned_target_sequence, realigned_offtarget_sequence


"""
Given a targetsite and window, use a fuzzy regex to align the targetsite to
the window. Returns the best match(es).
"""


def alignSequences(targetsite_sequence, window_sequence, edist=8,bulges=1,max_score=7):
    '''
    Taylor H update 12/08/2025 to take indels
    '''
    window_sequence = window_sequence.upper()
    query_regex = regexFromSequence(targetsite_sequence, indels = bulges ,errors=max_score)

    matches = list()

    for m in regex.finditer(query_regex, window_sequence, overlapped=True, flags=regex.REVERSE | regex.BESTMATCH):
        matches.append(("+",m))
    for m in regex.finditer(query_regex, reverseComplement(window_sequence), overlapped=True,flags=regex.REVERSE | regex.BESTMATCH):
        matches.append(("-", m))

    lowest_edit, lowest_mismatch = 100, 100
    chosen_alignment_m, chosen_alignment_strand_m = None, ''

    # Use regex to find the best match allowing only for mismatches
    for strand, alignment in matches:
        if alignment != None:
            mismatches, insertions, deletions = alignment.fuzzy_counts
            if deletions > 0:
                if alignment.span()[0] == alignment.fuzzy_changes[2][0] or alignment.span()[1] == alignment.fuzzy_changes[2][0]:
                    deletions=deletions-1
                    mismatches+=1
            edit = (mismatches+ insertions + deletions)
            if lowest_edit > edit and edist >= edit:  #store if lower edit distance
                chosen_alignment_m = alignment
                lowest_edit =edit
                chosen_alignment_strand_m = strand
                lowest_mismatch = mismatches
            elif lowest_edit == edit and edist >= edit:
                if mismatches > lowest_mismatch:#store if same edit distance & higher mismatches which means fewer bulges
                    chosen_alignment_m = alignment
                    lowest_edit = edit
                    chosen_alignment_strand_m = strand
                    lowest_mismatch = mismatches
            else:
                pass

    mismatches ,substitutions, insertions, deletions = 0,0,0,0
    chosen_alignment_strand_b = ''
    offtarget_sequence_no_bulge, start_no_bulge, end_no_bulge = '', '', ''
    bulged_offtarget_sequence, score, length, bulged_start,bulged_end, realigned_target = '', '', '', '', '', 'none'

    if chosen_alignment_m:
        mismatches, insertions, deletions = chosen_alignment_m.fuzzy_counts
        substitutions = mismatches
        chosen_alignment_strand_b = chosen_alignment_strand_m
        length = len(chosen_alignment_m.group())
        if insertions + deletions ==0:
            offtarget_sequence_no_bulge = chosen_alignment_m.group()
            start_no_bulge = chosen_alignment_m.start()
            end_no_bulge = chosen_alignment_m.end()
            chosen_alignment_strand_b = ''
        else:
            offtarget_sequence_no_bulge, start_no_bulge, end_no_bulge, chosen_alignment_strand_m = '','', '', ''


        if insertions + deletions >0:
            realigned_target, bulged_offtarget_sequence = realignedSequences(targetsite_sequence, chosen_alignment_m,
                                                                         max_score)
            score = substitutions + (insertions + deletions) * 3
            bulged_start = chosen_alignment_m.start()
            bulged_end = chosen_alignment_m.end()

    return [offtarget_sequence_no_bulge, mismatches, length, chosen_alignment_strand_m,
            start_no_bulge, end_no_bulge,
            realigned_target,
            bulged_offtarget_sequence, length, score, substitutions, insertions, deletions, chosen_alignment_strand_b,
            bulged_start, bulged_end]



""" Get sequences from some reference genome
"""


def get_sequence(reference_genome, chromosome, start, end, strand="+"):
    if strand == "+":
        seq = reference_genome[chromosome][int(start):int(end)]
    elif strand == "-":
        seq = reference_genome[chromosome][int(start):int(end)].reverse.complement
    return str(seq)


def bam_to_dict(bam, MAPQ=0):
    """Store all reads, paired/single, into dict
    """
    output = {}
    for read in bam.fetch(region=None):  # region='chr11:46408000-46409000'):
        qname = read.query_name
        if read.is_paired:
            if read.is_unmapped:
                continue
            if read.mapping_quality < MAPQ:
                continue
            if not qname in output:
                output[qname] = {0: [], 1: []}
            if read.is_read1:
                output[qname][0].append(read)
            else:
                output[qname][1].append(read)
    return output


def compare(ref, bam, control, targetsite, BEsearch_radius,BEmodel_min_overlap, BEmodel_max_overlap, window_size, mapq_threshold,
             mismatch_threshold, edist_threshold, bulge_threshold, name,
            out,  read_count_cutoff = 6):
    # housekeeping variables


    reference_genome_pyfaidx = pyfaidx.Fasta(ref)
    pkl_data = {'total nuclease count':None,
                'total control count': None,
                'total nuclease outie overlap':None,
                'total control outie overlap':None,
                'total nuclease matched':0,
                'total control matched': 0,
                'total nuclease unmatched':0,
                'total control unmatched': 0,
                }
    output_list = list()
    bg_position = list()  # List to store nuclease_position_counts that were observed at least once
    bg_narrow = list()  # List to store the sum of nuclease_position_counts in the narrow window

    output_dir = os.path.dirname(out)

    combined_ga = HTSeq.GenomicArray("auto", stranded=False)  # Store the union of control and nuclease positions
    offtarget_ga_windows = HTSeq.GenomicArray("auto", stranded=False)  # Store potential off-target sites
    ga_narrow_windows = HTSeq.GenomicArray("auto",
                                           stranded=False)  # Store potential off-target sites narrow windows read counts
    ga_narrow_windows_converted = HTSeq.GenomicArray("auto", stranded=False)
    ga_narrow_windows_noise = HTSeq.GenomicArray("auto", stranded=False)

    control_offtarget_ga_windows= HTSeq.GenomicArray("auto", stranded=False)
    control_ga_narrow_windows = HTSeq.GenomicArray("auto", stranded=False)
    control_ga_narrow_windows_converted = HTSeq.GenomicArray("auto", stranded=False)
    control_ga_narrow_windows_noise = HTSeq.GenomicArray("auto", stranded=False)
    #met critera,single end,base edited, noise (SE or PE)

    nuclease_ga, nuclease_ga_coverage_single, nuclease_ga_converted, nuclease_ga_noise, total_nuclease_count,total_nuclease_not_noise_count = \
        tabulate_start_positions_BE(bam=bam, label=name, output_dir=output_dir,mapq_threshold=mapq_threshold,
                                    BEmodel_min_overlap=BEmodel_min_overlap, BEmodel_max_overlap=BEmodel_max_overlap)

    control_ga, control_ga_coverage_single, control_ga_converted, control_ga_noise, total_control_count,total_control_not_noise_count = \
        tabulate_start_positions_BE(bam=control, label="Control_" + name, output_dir=output_dir,mapq_threshold=mapq_threshold,
                                    BEmodel_min_overlap=BEmodel_min_overlap, BEmodel_max_overlap=BEmodel_max_overlap)

    pkl_data['total nuclease count'] = total_nuclease_count
    pkl_data['total control count'] = total_control_count
    pkl_data['total nuclease outie overlap'] = total_nuclease_not_noise_count
    pkl_data['total control outie overlap'] = total_control_not_noise_count

    logger.info("Finished tabulating read start")
    # For all positions with detected read mapping positions, put into a combined genomicArray
    for iv, value in nuclease_ga.steps():
        if value:
            combined_ga[iv] = 1
    for iv, value in control_ga.steps():
        if value:
            combined_ga[iv] = 1

    logger.info("Finished combined_ga")
    for iv, value in combined_ga.steps():
        if value:
            for position in iv.range(step=1):
                # Define the windows

                window = HTSeq.GenomicInterval(position.chrom, max(0, position.pos - window_size),
                                               position.pos + window_size + 1)
                # Start mapping positions, at the specific base position
                nuclease_position_counts = nuclease_ga[position]

                control_position_counts = control_ga[position]
                # Store control_position_counts for which it was observed at least one read
                if control_position_counts >= 0:
                    bg_position.append(control_position_counts)

                # In the narrow (parameter-specified) window
                nuclease_window_counts = sum(nuclease_ga[window])
                control_window_counts = sum(control_ga[window])

                # new vars
                nuclease_window_converted_counts = sum(nuclease_ga_converted[window])
                nuclease_window_noise_counts = sum(nuclease_ga_noise[window])

                control_window_converted_counts = sum(control_ga_converted[window])
                control_window_noise_counts = sum(control_ga_noise[window])
                # print (position.chrom,position.pos,nuclease_position_counts,nuclease_window_converted_counts,nuclease_window_noise_counts)
                # Store control_window_counts greater than zero
                if control_window_counts >= 0:
                    bg_narrow.append(control_window_counts)

                # A list of the outputs
                row = [position.chrom, position.pos, nuclease_position_counts, control_position_counts,
                       nuclease_window_counts, control_window_counts, nuclease_window_converted_counts,
                       nuclease_window_noise_counts,control_window_converted_counts,
                       control_window_noise_counts]
                # nuclease_position_counts = count that are extactly this, nuclease_window_counts = counts in this wondow, nuclease_window_noise_counts did not meet critera
                output_list.append(row)


    for idx, fields in enumerate(output_list):

        if fields[2] >= read_count_cutoff or fields[4] >= read_count_cutoff:  # fields[2] is nuclease_position_counts and fields[4] is nuclease_window_counts, 4 should be alwasy higher
            read_chr = fields[0]
            read_position = fields[1]
            offtarget_ga_windows[HTSeq.GenomicPosition(read_chr, read_position, '.')] = 1
            ga_narrow_windows[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[
                4]  # this is the sum of number reads around a given position
            ga_narrow_windows_converted[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[
                6]  # this is the sum of number reads around a given position
            ga_narrow_windows_noise[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[
                7]  # this is the sum of number reads around a given position

            control_offtarget_ga_windows[HTSeq.GenomicPosition(read_chr, read_position, '.')] = 1
            control_ga_narrow_windows[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[
                5]  # this is the sum of number reads around a given positio
            control_ga_narrow_windows_converted[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[
                8]  # this is the sum of number reads around a given position
            control_ga_narrow_windows_noise[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[
                9]  # this is the sum of number reads around a given position
        # elif fields[4] >= edited_read_cutoff * 2 and fields[6] >= edited_read_cutoff: # this are base edited reads
        #     # print ("using edited reads")
        #     read_chr = fields[0]
        #     read_position = fields[1]
        #     offtarget_ga_windows[HTSeq.GenomicPosition(read_chr, read_position, '.')] = 1
        #     ga_narrow_windows[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[
        #         4]  # this is the sum of number reads around a given position
        #     ga_narrow_windows_converted[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[
        #         6]  # this is the sum of number reads around a given position
        #     ga_narrow_windows_noise[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[
        #         7]  # this is the sum of number reads around a given position


    ga_consolidated_windows = find_windows(offtarget_ga_windows, window_size)  # consolidate windows within 3 bp
    logger.info(f"Start off-target sequence fuzzy alignment using mismatch cutoff of {mismatch_threshold}")

    output_alignments(narrow_ga=ga_narrow_windows,
                      ga_windows=ga_consolidated_windows,
                      narrow_ga_converted = ga_narrow_windows_converted,
                      ga_narrow_windows_noise=ga_narrow_windows_noise,
                      control_narrow_ga=control_ga_narrow_windows,
                      control_narrow_ga_converted=control_ga_narrow_windows_converted,
                      control_ga_narrow_windows_noise=control_ga_narrow_windows_noise,
                      reference_genome=reference_genome_pyfaidx,
                      target_sequence=targetsite,
                      target_name=name,
                      bam_filename=bam,
                      edist_threshold=edist_threshold,
                      mismatch_threshold=mismatch_threshold,
                      bulge_threshold=bulge_threshold,
                      BEsearch_radius=BEsearch_radius,
                      pkl_data=pkl_data,
                      out= f"{output_dir}/{name}")

def main():
    parser = argparse.ArgumentParser(
        description='Identify off-target candidates from Illumina short read sequencing data.')
    parser.add_argument('--ref', help='Reference Genome Fasta', required=True)
    parser.add_argument('--bam', help='Sorted BAM file', required=True)
    parser.add_argument('--control', help='Control BAM file', required=True)
    parser.add_argument('--targetsite', help='Targetsite Sequence', required=True)
    parser.add_argument('--search_radius', help='Search radius around the position window', default=20, type=int)
    parser.add_argument('--window_size', help='Windowsize', default=3, type=int)
    parser.add_argument('--mapq', help='mapq threshold', default=50, type=int)
    parser.add_argument('--gap', help='Gap threshold', default=3, type=int)
    parser.add_argument('--start', help='Start threshold', default=1, type=int)
    parser.add_argument('--mismatch_threshold', help='Maximum score threshold', default=6, type=int)
    parser.add_argument('--read_count_cutoff', help='read_count threshold', default=6, type=int)
    parser.add_argument('--read_length', help='read_length', default=151, type=int)
    parser.add_argument('--merged', dest='merged', action='store_true', default=True)
    parser.add_argument('--all_chromosomes', dest='all_chromosomes', action='store_true', default=False)
    parser.add_argument('--name', help='Targetsite Name', required=False)
    parser.add_argument('--cells', help='Cells', required=False)
    parser.add_argument('--out', help='Output file base', required=True)
    args = parser.parse_args()


# Run the comparison if the control bam is specified, otherwise run the standard site identification routine.

# compare(args.ref, args.bam, args.control, args.targetsite, args.search_radius, args.window_size, args.mapq, args.gap,
# 		args.start, args.mismatch_threshold, args.name, args.cells, args.out, args.all_chromosomes, args.merged,args.read_count_cutoff)

if __name__ == "__main__":
    main()

