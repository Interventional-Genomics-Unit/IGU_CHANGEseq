from __future__ import print_function

import argparse
import HTSeq
import os
import pyfaidx
import regex
import pysam
import pandas as pd
from statsmodels.distributions.empirical_distribution import ECDF
import sys
import numpy as np
import logging
import pickle
""" Tabulate merged start positions.
	Identify genomic coordinates for reads mapping across 151/152 bp position.
	Add positions to genomic array.
"""
logger = logging.getLogger('root')
logger.propagate = False

def is_outie(x,y):
	if x.is_reverse:
		if not y.is_reverse:
			if x.reference_start < y.reference_start:
				return True,x.reference_end - y.reference_start
	else:
		if y.is_reverse:
			if y.reference_start < x.reference_start:
				return True,y.reference_end - x.reference_start
	return False,0



# def tabulate_merged_start_positions(bam,mapq_threshold, gap_threshold, label, start_threshold, output_dir, read_length=151):
# 	PE_read_stats_file = f"{output_dir}/{label}.PE.read_stats.csv"
# 	PE_read_stats_list_of_dict = []  # # read, chr, pos1, pos2, is_outie, overlap_bp, converted_pos_list
#
# 	logger.info("Reading bam to a dictionary")
# 	reads_dict = bam_to_dict_merged(pysam.AlignmentFile(bam, "rb", threads=48), mapq_threshold)
# 	n_reads = len(reads_dict.keys())
#
# 	logger.info("Finished bam to dictionary")
#
# 	### TH use bam writer to subset the bam file with reads that matches the criteria mapq
# 	#bam_writer = HTSeq.BAM_Writer.from_BAM_Reader(filtered_bam_filename,sorted_bam_file)
# 	ga_windows = HTSeq.GenomicArray("auto", stranded=False, typecode="d")
# 	ga_coverage = HTSeq.GenomicArray("auto", stranded=False, typecode="d")
# 	ga_noise = HTSeq.GenomicArray("auto", stranded=False, typecode="d")
# 	read_count = 0
# 	PE_read_count = 0
# 	not_noise_count = 0
#
# 	for read in reads_dict:
# 		noise_flag = True  # only count noise once
# 		read_count += 1
# 		if not read_count % 150000:
# 			logger.info(f'processing {read_count} ; {round(read_count / n_reads, 4) * 100}%')
# 		read1_list = reads_dict[read][0]
# 		read2_list = reads_dict[read][1]
#
# 		ga_coverage[read.iv] += 1
#
#
# 			for read in sorted_bam_file:
# 				output = False
# 				first_read_chr, first_read_position, first_read_strand = None, None, None
# 				second_read_chr, second_read_position, second_read_strand = None, None, None
# 				if read.aQual > mapq_threshold and read.aligned:
#
# 		for cigar_operation in read.cigar:
# 			# Identify positions that end in position 151 and start at position 151
# 			# Note strand polarity is reversed for position 151 (because it is part of the strand that was
# 			# reverse complemented initially before alignment
# 			# Yichao , replace 151 with read_length, replace 146, 156 as +-5, 5/12/2020
# 			if cigar_operation.type == 'M':
# 				if ((cigar_operation.query_from <= (read_length-5) - start_threshold) and
# 						(read_length - start_threshold <= cigar_operation.query_to)):
# 					first_read_cigar = cigar_operation
# 					first_read_chr = cigar_operation.ref_iv.chrom
# 					first_end = min(cigar_operation.query_to, read_length)
# 					distance = first_end - cigar_operation.query_from
# 					first_read_position = cigar_operation.ref_iv.start + distance - 1
# 					first_read_strand = '-'
# 				if ((cigar_operation.query_from <= read_length + start_threshold) and
# 						((read_length+5)  + start_threshold <= cigar_operation.query_to)):
# 					second_read_cigar = cigar_operation
# 					second_read_chr = cigar_operation.ref_iv.chrom
# 					second_end = max(read_length, cigar_operation.query_from)
# 					distance = second_end - cigar_operation.query_from
# 					second_read_position = cigar_operation.ref_iv.start + distance
# 					second_read_strand = '+'
#
# 		if first_read_chr:
# 			if first_read_chr == second_read_chr and (pattern.match(str(first_read_chr))) and \
# 							first_read_position is not None and second_read_position is not None:
# 				if abs(first_read_position - second_read_position) <= gap_threshold:
# 					output = True
# 					ga[HTSeq.GenomicPosition(first_read_chr, first_read_position, first_read_strand)] += 1
# 					ga_windows[HTSeq.GenomicPosition(first_read_chr, first_read_position, first_read_strand)] = 1
# 					#ga_stranded[HTSeq.GenomicPosition(first_read_chr, first_read_position, first_read_strand)] += 1
#
# 					ga[HTSeq.GenomicPosition(second_read_chr, second_read_position, second_read_strand)] += 1
# 					ga_windows[HTSeq.GenomicPosition(second_read_chr, second_read_position, second_read_strand)] = 1
# 					#ga_stranded[HTSeq.GenomicPosition(second_read_chr, second_read_position, second_read_strand)] += 1
# 					#bam_writer.write(read)
#
# 		if output == True:
# 			print(name, targetsite, cells, filename_base, first_read_chr, first_read_position,
# 				  first_read_strand, second_read_chr, second_read_position, second_read_strand, sep='\t', file=o)
#
# 	read_count += 1
#
#
# 	return ga, ga_windows, ga_coverage, read_count


""" Tabulate the start positions for the 2nd read in pair across the genome.
	Only consider alignments with matching positions from the beginning of the read.
	For read pairs with multiple alignments, pick the one with matching positions at the beginning.
"""
# def bam_to_dict_merged(bam, MAPQ):
# 	"""Store all reads, paired/single, into dict
# 	"""
# 	output = {}
# 	for read in bam.fetch(region=None):  # region='chr11:46408000-46409000'):
# 		qname = read.query_name
#
# 		if read.is_paired:
# 			if read.is_unmapped:
# 				continue
# 			if read.mapping_quality < MAPQ:
# 				continue
# 			if not qname in output:
# 				output[qname] = {0: [], 1: []}
# 			if read.is_read1:
# 				output[qname][0].append(read)
# 			else:
# 				output[qname][1].append(read)
# 	return output
def bam_to_dict(bam, MAPQ):
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

def get_read_strand(read):
	if read.is_reverse:
		return "-"
	return "+"


def get_read_start(read):
	if read.is_reverse:
		return read.reference_end
	return read.reference_start + 1


def tabulate_start_positions(bam,mapq_threshold,
							 gap_threshold, label,output_dir):
	PE_read_stats_file = f"{output_dir}/{label}.PE.read_stats.csv"
	PE_read_stats_list_of_dict = []  # # read, chr, pos1, pos2, is_outie, overlap_bp, converted_pos_list

	logger.info("Reading bam to a dictionary")
	reads_dict = bam_to_dict(pysam.AlignmentFile(bam, "rb", threads=48), mapq_threshold)
	n_reads = len(reads_dict.keys())

	logger.info("Finished bam to dictionary")

	ga_windows = HTSeq.GenomicArray("auto", stranded=False,typecode="d")
	ga= HTSeq.GenomicArray("auto", stranded=False,typecode="d")
	ga_noise = HTSeq.GenomicArray("auto", stranded=False, typecode="d")
	read_count = 0
	PE_read_count = 0
	not_noise_count = 0

	for read in reads_dict:

		read_count += 1
		if not read_count % 150000:
			logger.info(f'processing {read_count} ; {round(read_count / n_reads, 4) * 100}%')
		read1_list = reads_dict[read][0]
		read2_list = reads_dict[read][1]


		if len(read1_list) == 0 or len(read2_list) == 0:
			pass

		else:
			PE_read_count += 1
			for i in read1_list:
				read1_start = get_read_start(i)
				read1_strand = get_read_strand(i)
				for j in read2_list:
					read2_start = get_read_start(j)
					read2_strand = get_read_strand(j)

					if i.reference_name == j.reference_name:
						flag, overlap_bp = is_outie(i, j)
						if flag:
							noise_flag = True
							if abs(overlap_bp) >= 0 and abs(overlap_bp) <= gap_threshold:
								noise_flag = False
						else:
							noise_flag = True


			if noise_flag:
				ga_noise[HTSeq.GenomicPosition(i.reference_name, i.reference_start)] += 1
				ga_noise[HTSeq.GenomicPosition(j.reference_name, j.reference_start)] += 1
			else:

				ga[HTSeq.GenomicPosition(i.reference_name, read1_start)] += 1
				ga[HTSeq.GenomicPosition(i.reference_name, read2_start)] += 1
				ga_windows[HTSeq.GenomicPosition(i.reference_name, read1_start)] = 1
				ga_windows[HTSeq.GenomicPosition(i.reference_name, read2_start)] = 1
				not_noise_count += 1
				PE_read_stats_list_of_dict.append({"read": read,
												   "chr": i.reference_name,
												   "pos1": i.reference_start,
												   "strand1": read1_strand,
												   "pos2": j.reference_start,
												   "strand2": read2_strand,
												   "is_outie": flag,
												   "overlap_bp": overlap_bp})
	pd.DataFrame.from_dict(PE_read_stats_list_of_dict).to_csv(PE_read_stats_file, index=False)

	logger.info(f"{len(reads_dict.keys())} reads were mapped and met MAPQ threshold: {mapq_threshold}")
	logger.info(f"{not_noise_count}, {round(not_noise_count / n_reads, 4) * 100}% outie paired reads <= {gap_threshold}")

	return ga ,ga_windows, ga_noise, PE_read_count,not_noise_count



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

def output_alignments(narrow_ga, ga_windows, ga_narrow_windows_noise,control_narrow_ga,
					  control_ga_narrow_windows_noise,reference_genome,
					  target_sequence,
					  target_name,
					  bam_filename,
					  edist_threshold,mismatch_threshold,bulge_threshold, search_radius, pkl_data,out):
	pkl_out = '{0}_total_counts.pkl'.format(out)
	outfile_matched = '{0}_identified_matched.txt'.format(out)
	outfile_unmatched = '{0}_identified_unmatched.txt'.format(out)

	# dictionary to store the matched reads
	matched_dict = {}
	# dictionary to add read_count for each pair chromosome:start_position among matched reads
	reads_dict = {}
	# dictionary to store window_start. For duplicated matched off-target.
	control_reads_dict = {}
	window_min = {}
	# dictionary to store window_end. For duplicated matched off-target.
	window_max = {}

	# dictionary to store the unmatched reads
	unmatched_dict = {}

	for iv, value in ga_windows.steps():
		if value:
			window_sequence = get_sequence(reference_genome, iv.chrom, iv.start - search_radius, iv.end + search_radius)

			(start, end, strand, offtarget_sequence, offtarget_sequence_aligned, target_aligned, length, mismatches, insertions, deletions,
			 alignment_found,bulge_alignment_flag) = \
				alignSequences(target_sequence, window_sequence, edist=edist_threshold, bulges=bulge_threshold,
							   max_score=mismatch_threshold)

			# get genomic coordinates of sequences
			mm_start, mm_end, b_start, b_end = '', '', '', ''
			if alignment_found:
				if strand == '+':
					mm_start = iv.start - search_radius + int(start)
					mm_end = iv.start - search_radius + int(end)
				if strand == '-':
					mm_start = iv.end + search_radius - int(end)
					mm_end = iv.end + search_radius - int(start)


			#  define overall start, end and strand. For bed annotation purposes
			if alignment_found:
				target_start_absolute = mm_start
				target_end_absolute = mm_end
				target_strand_absolute = strand
			else:
				target_start_absolute = iv.start
				target_end_absolute = iv.end
				target_strand_absolute = '*'

			name = iv.chrom + ':' + str(target_start_absolute) + '-' + str(target_end_absolute)
			read_count = int(max(set(narrow_ga[iv])))
			control_read_count = int(max(set(control_narrow_ga[iv])))
			nuclease_noise_count = str(int(max(ga_narrow_windows_noise[iv])))
			control_noise_count = str(int(max(control_ga_narrow_windows_noise[iv])))
			filename = os.path.basename(bam_filename)
			full_name = str(target_name) + '_' + str(name)

			if alignment_found:
				tag = iv.chrom + ':' + str(target_start_absolute)
				if tag not in reads_dict.keys():
					reads_dict[tag] = read_count
					control_reads_dict[tag] = control_read_count
					window_min[tag] = [iv.start]
					window_max[tag] = [iv.end]
					matched_dict[tag] = [iv.chrom, target_start_absolute, target_end_absolute, name,target_strand_absolute,
										 read_count,control_read_count,
										 offtarget_sequence, offtarget_sequence_aligned,
										 mismatches, insertions, deletions,
										 filename, target_name, target_name, full_name, target_sequence,
										 target_aligned,0, 0,nuclease_noise_count,control_noise_count, iv.start, iv.end, iv, window_sequence]

				else:
					current_read_count = reads_dict[tag]
					current_control_read_count = control_reads_dict[tag]
					reads_dict[tag] = max(current_read_count, read_count)
					control_reads_dict[tag] = max(current_control_read_count, control_read_count)
					window_min[tag].append(iv.start)
					window_max[tag].append(iv.end)
					matched_dict[tag] = [iv.chrom, target_start_absolute, target_end_absolute, name,target_strand_absolute,
										 reads_dict[tag], control_reads_dict[tag],
										 offtarget_sequence, offtarget_sequence_aligned,
										 mismatches, insertions, deletions,
										 filename, target_name, target_name, full_name, target_sequence,
										 target_aligned, 0, 0,
										 nuclease_noise_count, control_noise_count,min(window_min[tag]), max(window_max[tag]), iv, window_sequence]

			# no matching site but meets other critera
			else:
				untag = iv.chrom + ':' + str(iv.start)
				unmatched_dict[untag] = [iv.chrom, target_start_absolute, target_end_absolute, name,
										 target_strand_absolute,
										 read_count, control_read_count,
										 offtarget_sequence, offtarget_sequence_aligned,
										 mismatches, insertions, deletions,
										 filename, target_name, target_name, full_name, target_sequence,
										 target_aligned, 0, 0,
										 nuclease_noise_count, control_noise_count]


	# Yichao, add control reads
	logger.info("Writing matched table")
	tags_sorted = matched_dict.keys()
	tags_sorted = sorted(tags_sorted)

	pkl_data['total nuclease matched'] = [x for x in matched_dict.values()]
	pkl_data['total nuclease matched'] = len(tags_sorted)

	o1 = open(outfile_matched, 'w')

	print('#Chromosome', 'Start', 'End', 'Genomic Coordinate', 'Strand',
		  'Nuclease_Read_Count',  'Control_Read_Count',
		  'Site_Sequence', 'Site_Sequence_Gaps_Allowed',
		  'Site_Substitution_Number', 'RNA_Bulge','DNA_Bulge',
		  'File_Name', 'Cell', 'Target_site', 'Full_Name', 'Target_Sequence', 'Realigned_Target_Sequence',
		  'Nuclease_Base_Converted', 'Control_Base_Converted', 'Nuclease_Noise', 'Control_Noise',
		  'MappingPositionStart', 'MappingPositionEnd',
		  'WindowName', 'WindowSequence',  # which column, -2 -1
		  sep='\t', file=o1)
	o1.close()

	with open(outfile_matched, 'a') as o1:
		for key in tags_sorted:
			row = matched_dict[key]
			pkl_data['total nuclease matched'] += int(row[5])
			pkl_data['total control matched'] += int(row[6])
			print(*(row), sep='\t', file=o1)

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
	bulge_alignment_flag = False
	alignment_found = False

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

	mismatches ,insertions, deletions = 0,0,0
	offtarget_sequence, start, end,strand = '', '', '',''
	offtarget_sequence_aligned, length, target_aligned = '', '', '',

	if chosen_alignment_m:
		alignment_found = True
		mismatches, insertions, deletions = chosen_alignment_m.fuzzy_counts
		length = len(chosen_alignment_m.group())
		strand = chosen_alignment_strand_m
		offtarget_sequence = chosen_alignment_m.group()
		start = chosen_alignment_m.start()
		end = chosen_alignment_m.end()
		if insertions + deletions !=0:
			bulge_alignment_flag = True
			target_aligned,offtarget_sequence_aligned = realignedSequences(targetsite_sequence, chosen_alignment_m,
																		 max_score)
		else:
			target_aligned, offtarget_sequence_aligned = offtarget_sequence, targetsite_sequence


	#return [offtarget_sequence_no_bulge, mismatches, length, chosen_alignment_strand_m,
	#		start_no_bulge, end_no_bulge,
	#		realigned_target,
	#		bulged_offtarget_sequence, length, score, substitutions, insertions, deletions, chosen_alignment_strand_b,
	#		bulged_start, bulged_end]
	return [start,end,strand,offtarget_sequence, offtarget_sequence_aligned,target_aligned, length,mismatches, insertions, deletions,alignment_found,bulge_alignment_flag]

""" Get sequences from some reference genome
"""
def get_sequence(reference_genome, chromosome, start, end, strand="+"):
	if strand == "+":
		seq = reference_genome[chromosome][int(start):int(end)]
	elif strand == "-":
		seq = reference_genome[chromosome][int(start):int(end)].reverse.complement
	return str(seq)


def compare(ref, bam, control, targetsite, search_radius, windowsize, mapq_threshold, gap_threshold,
			start_threshold, mismatch_threshold,edist_threshold,bulge_threshold, name,
			out,  merged=False,read_count_cutoff=6,read_length=151):


	reference_genome = pyfaidx.Fasta(ref)

	pkl_data = {'total nuclease count': None,
				'total control count': None,
				'total nuclease outie overlap': None,
				'total control outie overlap': None,
				'total nuclease matched': 0,
				'total control matched': 0,
				'total nuclease unmatched': 0,
				'total control unmatched': 0,
				}

	combined_ga = HTSeq.GenomicArray("auto", stranded=False)  # Store the union of control and nuclease positions
	offtarget_ga_windows = HTSeq.GenomicArray("auto", stranded=False)  # Store potential off-target sites
	ga_narrow_windows = HTSeq.GenomicArray("auto", stranded=False)  # Store potential off-target sites narrow windows read counts
	ga_narrow_windows_noise = HTSeq.GenomicArray("auto", stranded=False)

	control_offtarget_ga_windows = HTSeq.GenomicArray("auto", stranded=False)
	control_ga_narrow_windows = HTSeq.GenomicArray("auto", stranded=False)
	control_ga_narrow_windows_noise = HTSeq.GenomicArray("auto", stranded=False)

	
	output_list = list()
	bg_position = list()  # List to store nuclease_position_counts that were observed at least once
	bg_narrow = list()  # List to store the sum of nuclease_position_counts in the narrow window

	output_dir = os.path.dirname(out)


	logger.info("Tabulate nuclease standard start positions.")

	nuclease_ga, nuclease_ga_windows, nuclease_ga_noise, total_nuclease_count,total_nuclease_not_noise_count = \
		tabulate_start_positions(bam=bam,mapq_threshold = mapq_threshold,
								 gap_threshold = gap_threshold,label = name, output_dir=output_dir)
	logger.info("Tabulate control standard start positions.")
	control_ga, control_ga_windows, control_ga_noise, total_control_count,total_control_not_noise_count = \
		tabulate_start_positions(bam=control, mapq_threshold =mapq_threshold,
								 gap_threshold =  gap_threshold, label = "Control_" + name, output_dir=output_dir)

	pkl_data['total nuclease count'] = total_nuclease_count
	pkl_data['total control count'] = total_control_count
	pkl_data['total nuclease outie overlap'] = total_nuclease_not_noise_count
	pkl_data['total control outie overlap'] = total_control_not_noise_count

	output_filename = out + '_count.txt'
	with open(output_filename, 'w') as o:
		logger.info("Writing counts to {0}".format(output_filename))

		# For all positions with detected read mapping positions, put into a combined genomicArray
		for iv, value in nuclease_ga.steps():
			if value:
				combined_ga[iv] = 1
		for iv, value in control_ga.steps():
			if value:
				combined_ga[iv] = 1

		for iv, value in combined_ga.steps():
			if value:
				for position in iv.range(step=1):


					# Define the windows
					window = HTSeq.GenomicInterval(position.chrom, max(0, position.pos - windowsize),
												   position.pos + windowsize + 1)


					nuclease_position_counts = nuclease_ga[position]
					control_position_counts = control_ga[position]
					# Store control_position_counts for which it was observed at least one read
					if control_position_counts >= 0:
						bg_position.append(control_position_counts)

					# In the narrow (parameter-specified) window
					nuclease_window_counts = max(nuclease_ga[window])
					control_window_counts = sum(control_ga[window])

					nuclease_window_noise_counts = sum(nuclease_ga_noise[window])
					control_window_noise_counts = sum(control_ga_noise[window])

					# Store control_window_counts greater than zero
					if control_window_counts >= 0:
						bg_narrow.append(control_window_counts)

					# A list of the outputs
					row = [position.chrom, position.pos, nuclease_position_counts, control_position_counts,
						   nuclease_window_counts, control_window_counts,nuclease_window_noise_counts,control_window_noise_counts]
					output_list.append(row)

		print('#Chromosome', 'zero_based_Position', 'Nuclease_Position_Reads', 'Control_Position_Reads',
			  'Nuclease_Window_Reads', 'Control_Window_Reads', 'Nuclease_Noise_Reads', 'Control_Noise_Reads', file=o, sep='\t')

		for idx, fields in enumerate(output_list):
			## read threshold is here: Yichao, read_count_cutoff
			if fields[2] >= read_count_cutoff or fields[4] >= read_count_cutoff:  # fields[2] is nuclease_position_counts and fields[4] is nuclease_window_counts
				read_chr = fields[0]
				read_position = fields[1]
				offtarget_ga_windows[HTSeq.GenomicPosition(read_chr, read_position, '.')] = 1
				ga_narrow_windows[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[4]
				ga_narrow_windows_noise[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[
					6]  # this is the sum of number reads around a given position


				control_offtarget_ga_windows[HTSeq.GenomicPosition(read_chr, read_position, '.')] = 1
				control_ga_narrow_windows[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[
					5]  # this is the sum of number reads around a given positio
				control_ga_narrow_windows_noise[HTSeq.GenomicPosition(read_chr, read_position, '.')] = fields[
					7]  # this is the sum of number reads around a given position


			print(*(fields),file=o, sep='\t')


		ga_consolidated_windows = find_windows(offtarget_ga_windows, windowsize)    # consolidate windows within 3 bp


		output_alignments(narrow_ga=ga_narrow_windows,
						  ga_windows = ga_consolidated_windows,
						  ga_narrow_windows_noise = ga_narrow_windows_noise,
						  control_narrow_ga=control_ga_narrow_windows,
						  control_ga_narrow_windows_noise=control_ga_narrow_windows_noise,
						  reference_genome=reference_genome,
						  target_sequence=targetsite,
						  target_name=name,
						  bam_filename=bam,
						  edist_threshold=edist_threshold,
						  mismatch_threshold=mismatch_threshold,
						  bulge_threshold=bulge_threshold,
						  search_radius=search_radius,
						  pkl_data=pkl_data,
						  out=f"{output_dir}/{name}")


def main():
	parser = argparse.ArgumentParser(description='Identify off-target candidates from Illumina short read sequencing data.')
	parser.add_argument('--ref', help='Reference Genome Fasta', required=True)
	parser.add_argument('--bam', help='Sorted BAM file', required=True)
	parser.add_argument('--control', help='Control BAM file', required=True)
	parser.add_argument('--targetsite', help='Targetsite Sequence', required=True)
	parser.add_argument('--search_radius', help='Search radius around the position window', default=20, type=int)
	parser.add_argument('--windowsize', help='Windowsize', default=3, type=int)
	parser.add_argument('--mapq', help='mapq threshold', default=50, type=int)
	parser.add_argument('--gap', help='Gap threshold', default=3, type=int)
	parser.add_argument('--start', help='Start threshold', default=1 , type=int)
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
	print("Nuclease: {0}\nControl: {1}".format(args.bam, args.control), file=sys.stderr)
	compare(args.ref, args.bam, args.control, args.targetsite, args.search_radius, args.windowsize, args.mapq, args.gap,
			args.start, args.mismatch_threshold, args.name, args.cells, args.out, args.all_chromosomes, args.merged,args.read_count_cutoff)

if __name__ == "__main__":
	main()