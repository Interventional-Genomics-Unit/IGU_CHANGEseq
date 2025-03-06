from __future__ import print_function
import svgwrite
import os
import logging
import argparse
from Bio import SeqUtils


### 2017-October-11: Adapt plots to new output; inputs are managed using "argparse".

logger = logging.getLogger('root')
logger.propagate = False

boxWidth = 10
box_size = 15
v_spacing = 3

global colors
colors = ['#EC5121', '#96C04B', '#FFAC0F', '#7ACBC8', '#97FFFF']
#colors = {'G': '#00CDCD', 'A': '#3A5FCD', 'T': '#8EE5EE', 'C': '#B0E0E6', 'N': '#97FFFF',
#          'R': '#B3B3B3', '-': '#B3B3B3'}


def levenshtein_two_matrix_rows(str1, str2, PAM):
    # Get the lengths of the input strings
    str1 = str1[:-len(PAM)]
    str2 = str2[:-len(PAM)]

    dist = 0
    for i,bp in enumerate(str1):
        if str2[i] != bp:
            dist+=1

    return dist

def parseSitesFile(infile):
    offtargets = []
    total_seq = 0
    with open(infile, 'r') as f:
        f.readline()
        for line in f:
            line = line.rstrip('\n')
            if '\t' in line:
                line_items = line.split('\t')
            else:
                line_items = line.split(',')

            offtarget_reads = line_items[4]
            no_bulge_offtarget_sequence = line_items[7].upper()  # Site_Sequence
            bulge_offtarget_sequence = line_items[9]
            target_seq = line_items[16]
            realigned_target_seq = line_items[17]
            coord = line_items[3].split("-")[0]
            num_mismatch = int(line_items[8])

            if "intergenic" not in line_items[-1]:
                annot = line_items[-2] + "," + line_items[-1].replace("non-coding RNA","ncRNA") # gene name and feature
            else:
                annot = ""


            if no_bulge_offtarget_sequence != '' or bulge_offtarget_sequence != '':
                #print(no_bulge_offtarget_sequence,bulge_offtarget_sequence)
                if no_bulge_offtarget_sequence:
                    total_seq += 1
                if bulge_offtarget_sequence:
                    total_seq += 1
                    #num_mismatch = 0
                    #length = len(realigned_target_seq.strip())
                    #for l in range(length):
                    #    k = bulge_offtarget_sequence.strip()[l]
                    #    j = realigned_target_seq.strip()[l]
                    #    if k != '-' and j != '-' and k !=j:
                    #        num_mismatch += 1

                offtargets.append({'seq': no_bulge_offtarget_sequence.strip(),
                                   'bulged_seq': bulge_offtarget_sequence.strip(),
                                   'reads': int(offtarget_reads.strip()),
                                   'coord': str(coord),
                                   'annot': str(annot),
                                   'num_mismatch': str(num_mismatch),
                                   'target_seq': target_seq.strip(),
                                   'realigned_target_seq': realigned_target_seq.strip()
                                   })
    offtargets = sorted(offtargets, key=lambda x: x['reads'], reverse=True)

    return offtargets, target_seq, total_seq

# 3/6/2020 Yichao
def check_mismatch(a,b):
    from Bio.Data import IUPACData
    dna_dict = IUPACData.ambiguous_dna_values
    set_a = dna_dict[a.upper()]
    set_b = dna_dict[b.upper()]
    overlap = list(set(list(set_a)).intersection(list(set_b)))
    if len(overlap) == 0:
        return True
    else:
        return False

def find_PAM(seq,PAM):
    try:
        PAM_index = seq.index(PAM)
    except:# PAM on the left
        left_search = SeqUtils.nt_search(seq[:len(PAM)], PAM)
        if len(left_search)>1:
            PAM_index = left_search[1]
        else:
            right_search = SeqUtils.nt_search(seq[-len(PAM):], PAM)
            if len(right_search)>1:
                PAM_index = len(seq)-len(PAM)
            else:
                print ("PAM: %s not found in %s. Set PAM index to 20"%(PAM,seq))
                PAM_index=20
    return PAM_index
def draw_plot(target_seq,offtargets,total_seq,outfile,title,PAM):
    colors = {'G': '#EC5121', 'A': '#96C04B', 'T': '#FFAC0F', 'C': '#7ACBC8', 'N': '#97FFFF',
              'R': '#B3B3B3', '-': '#B3B3B3'}

    for c in ['Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', '.']:
        colors[c] = "#B3B3B3"
    # Initiate canvas
    dwg = svgwrite.Drawing(outfile, profile='full', size=(u'100%', 100 + total_seq * (box_size + 1)))

    if title is not None:
        # Define top and left margins
        x_offset = 20
        y_offset = 50
        dwg.add(dwg.text(title, insert=(x_offset, 30), style="font-size:20px; font-family:Courier"))
    else:
        # Define top and left margins
        x_offset = 20
        y_offset = 20

    tick_locations = []
    tick_legend = []
    # PAM_index = target_seq.index(PAM)
    PAM_index = find_PAM(target_seq, PAM)
    count = 0
    for i in range(PAM_index, 0, -1):
        count = count + 1
        if count % 10 == 0:
            tick_legend.append(count)
            # print (count,i)
            tick_locations.append(i)
    if len(PAM) >= 3:
        tick_legend += ['P', 'A', 'M'] + ['-'] * (len(PAM) - 3)
    else:
        tick_legend += ["PAM"] + ['-'] * (len(PAM) - 3)
    tick_locations += range(PAM_index + 1, len(target_seq) + 1)
    if PAM_index == 0:
        tick_legend = []
        tick_locations = []
        tick_legend += ['P', 'A', 'M'] + ['-'] * (len(PAM) - 3)
        tick_locations += range(1, len(PAM) + 1)
        count = 0
        for i in range(len(PAM) + 1, len(target_seq) + 1):
            count = count + 1
            if count % 10 == 0 or count == 1:
                tick_legend.append(count)
                # print (count,i)
                tick_locations.append(i)
    # print (zip(tick_locations, tick_legend))
    for x, y in zip(tick_locations, tick_legend):
        dwg.add(dwg.text(y, insert=(x_offset + (x - 1) * box_size + 2, y_offset - 2),
                         style="font-size:10px; font-family:Courier"))

    # Draw reference sequence row
    for i, c in enumerate(target_seq):
        y = y_offset
        x = x_offset + i * box_size
        dwg.add(dwg.rect((x, y), (box_size, box_size), opacity=0.5, fill=colors[c]))
        dwg.add(
            dwg.text(c, insert=(x + 3, y + box_size - 3), fill='black', style="font-size:15px; font-family:Courier"))
    dwg.add(dwg.text('Reads', insert=(x_offset + box_size * len(target_seq) + 16, y_offset + box_size - 3),
                     style="font-size:15px; font-family:Courier"))
    dwg.add(dwg.text('Distance', insert=(box_size * (len(target_seq) + 1) + 90, y_offset + box_size - 3),
                     style="font-size:15px; font-family:Courier"))
    dwg.add(dwg.text('Annotation', insert=(box_size * (len(target_seq) + 1) + 200, y_offset + box_size - 3),
                     style="font-size:15px; font-family:Courier"))
    dwg.add(dwg.text('Coordinates', insert=(box_size * (len(target_seq) + 1) + 390, y_offset + box_size - 3),
                     style="font-size:15px; font-family:Courier"))

    # Draw aligned sequence rows
    y_offset += 1  # leave some extra space after the reference row
    line_number = 0  # keep track of plotted sequences
    for j, seq in enumerate(offtargets):
        realigned_target_seq = offtargets[j]['realigned_target_seq']
        no_bulge_offtarget_sequence = offtargets[j]['seq']
        bulge_offtarget_sequence = offtargets[j]['bulged_seq']
        b_dist = ''

        if no_bulge_offtarget_sequence != '': ## Theres at least an alignment with nok bulges
            b_dist = levenshtein_two_matrix_rows(no_bulge_offtarget_sequence, target_seq,PAM)
            k = 0
            line_number += 1
            y = y_offset + line_number * box_size
            for i, (c, r) in enumerate(zip(no_bulge_offtarget_sequence, target_seq)):
                x = x_offset + k * box_size
                if r == '-':
                    if 0 < k < len(target_seq):
                        x = x_offset + (k - 0.25) * box_size
                        dwg.add(dwg.rect((x, box_size * 1.4 + y), (box_size * 0.6, box_size * 0.6), opacity=0.6,
                                         fill=colors[c]))
                        dwg.add(dwg.text(c, insert=(x + 1, 2 * box_size + y - 2), fill='black',
                                         style="font-size:10px; font-family:Courier"))
                elif c == r:
                    dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black',
                                     style="font-size:10px; font-family:Courier"))
                    k += 1
                elif r == 'N':
                    dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black',
                                     style="font-size:15px; font-family:Courier"))
                    k += 1
                else:
                    dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), opacity=0.6, fill=colors[c]))
                    dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black',
                                     style="font-size:15px; font-family:Courier"))
                    k += 1
        if bulge_offtarget_sequence != '': #theres a bulge seq present
            k = 0
            line_number += 1
            y = y_offset + line_number * box_size
            b_dist = levenshtein_two_matrix_rows(bulge_offtarget_sequence, realigned_target_seq,PAM)
            for i, (c, r) in enumerate(zip(bulge_offtarget_sequence, realigned_target_seq)):
                x = x_offset + k * box_size
                if r == '-':
                    if 0 < k < len(realigned_target_seq):
                        x = x_offset + (k - 0.25) * box_size
                        dwg.add(dwg.rect((x, box_size * 1.4 + y), (box_size * 0.6, box_size * 0.6), opacity=0.6,
                                         fill=colors[c]))
                        dwg.add(dwg.text(c, insert=(x + 1, 2 * box_size + y - 2), fill='black',
                                         style="font-size:10px; font-family:Courier"))
                elif c == r:
                    dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black',
                                     style="font-size:10px; font-family:Courier"))
                    k += 1
                elif r == 'N':
                    dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black',
                                     style="font-size:15px; font-family:Courier"))
                    k += 1
                else:
                    dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), opacity=0.6, fill=colors[c]))
                    dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black',
                                     style="font-size:15px; font-family:Courier"))
                    k += 1

        if no_bulge_offtarget_sequence == '' or bulge_offtarget_sequence == '': # there are not two different alignments
            reads_text = dwg.text(str(seq['reads']), insert=(
                box_size * (len(target_seq) + 1) + 20, y_offset + box_size * (line_number + 2) - 2),
                                  fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(reads_text)
            mismatch_text = dwg.text(b_dist, insert=(
                box_size * (len(target_seq) + 1) + 130, y_offset + box_size * (line_number + 2) - 2),
                                     fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(mismatch_text)

            annot_text = dwg.text(seq['annot'], insert=(
                box_size * (len(target_seq) + 1) + 200, y_offset + box_size * (line_number + 2) - 2),
                                  fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(annot_text)

            coord_text = dwg.text(seq['coord'], insert=(
                box_size * (len(target_seq) + 1) + 390, y_offset + box_size * (line_number + 2) - 2),
                                  fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(coord_text)

        else:
            reads_text = dwg.text(str(seq['reads']), insert=(
            box_size * (len(target_seq) + 1) + 20, y_offset + box_size * (line_number + 1) + 5),
                                  fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(reads_text)

            mismatch_text = dwg.text(b_dist, insert=(
                box_size * (len(target_seq) + 1) + 130, y_offset + box_size * (line_number + 2) - 2),
                                     fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(mismatch_text)
            annot_text = dwg.text(seq['annot'], insert=(
                box_size * (len(target_seq) + 1) + 200, y_offset + box_size * (line_number + 1) + 5),
                                  fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(annot_text)
            mismatch_text = dwg.text(seq['coord'], insert=(
                box_size * (len(target_seq) + 1) + 380, y_offset + box_size * (line_number + 1) + 5),
                                     fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(mismatch_text)

            dist_text = dwg.text(b_dist, insert=(
                box_size * (len(target_seq) + 1) + 130, y_offset + box_size * (line_number + 1) - 2),
                                     fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(dist_text)
            ## TH added coords

            reads_text02 = dwg.text(u"\u007D", insert=(
            box_size * (len(target_seq) + 1) + 7, y_offset + box_size * (line_number + 1) + 5),
                                    fill='black', style="font-size:23px; font-family:Courier")
            dwg.add(reads_text02)
    dwg.save()

def visualizeOfftargets(infile, outfile, title, PAM):

    output_folder = os.path.dirname(outfile)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Get offtargets array from file
    offtargets, target_seq, total_seq = parseSitesFile(infile)

    draw_plot(target_seq, offtargets, total_seq, outfile, title, PAM)





def main():
    parser = argparse.ArgumentParser(description='Plot visualization plots for re-aligned reads.')
    parser.add_argument("-f", "--identified_file", help="FullPath/output file from reAlignment_circleseq.py",
                        required=True)
    parser.add_argument("-o", "--outfile", help="FullPath/VIZ", required=True)
    parser.add_argument("-t", "--title", help="Plot title", required=True)
    parser.add_argument("--PAM", help="PAM sequence", default="NGG")
    args = parser.parse_args()

    print(args)

    visualizeOfftargets(args.identified_file, args.outfile, args.title, args.PAM)


if __name__ == "__main__":
    main()


'''
infile = "/groups/clinical/projects/Assay_Dev/CHANGEseq/CS_12/identified/py3_781_12878_rep1_identified_matched_annotated.csv"
outfile = "/groups/clinical/projects/Assay_Dev/CHANGEseq/CS_12/visualization/py3_781_12878_rep1."
title = 'py3_781_12878_rep1'
visualizeOfftargets(infile = infile, outfile = outfile, title = title, PAM= 'NGG')
'''
