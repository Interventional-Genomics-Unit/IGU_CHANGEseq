from __future__ import print_function
import svgwrite
import os
import logging
import argparse
import pandas as pd
from Bio import SeqUtils


### 2017-October-11: Adapt plots to new output; inputs are managed using "argparse".
logger = logging.getLogger('root')
logger.propagate = False

boxWidth = 10
box_size = 15
v_spacing = 3

global colors
colors = ['#EC5121', '#96C04B', '#FFAC0F', '#7ACBC8', '#97FFFF']

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

    if infile.endswith(".txt"):
        df = pd.read_csv(infile, sep='\t')
    else:
        df = pd.read_csv(infile)
    for i,row in df.iterrows():
        bulge_flag = True

        offtarget_reads = int(row['Nuclease_Read_Count'])
        no_bulge_offtarget_sequence = row['Site_Sequence'].upper()  # Site_Sequence
        target_seq = row['Target_Sequence'].replace("-","")
        realigned_target_seq = row['Target_site']
        coord = row['Genomic Coordinate']
        num_mismatch = str(row['Site_Substitution_Number'])
        dist = int(row['RNA_Bulge']) + int(row['RNA_Bulge']) + int(row['Site_Substitution_Number'])

        if "intergenic" not in row['Feature']:
            annot =  row['Gene_Name'] + "," + row.Feature.replace("non-coding RNA","ncRNA") # gene name and feature
        else:
            annot = ""

        total_seq += 1

        if row['RNA_Bulge'] + row['RNA_Bulge'] == 0:
            bulge_flag = False

        offtargets.append({'seq': no_bulge_offtarget_sequence.strip(),
                           'reads': offtarget_reads,
                           'dist' : dist,
                           'coord': str(coord),
                           'annot': str(annot),
                           'num_mismatch': num_mismatch,
                           'target_seq': target_seq.strip(),
                           'realigned_target_seq': realigned_target_seq.strip(),
                           'bulge_flag': bulge_flag})
    offtargets = sorted(offtargets, key=lambda x: x['reads'], reverse=True)

    return offtargets, target_seq, total_seq


def find_PAM(seq,PAM):
    try:
        PAM_index = seq.index(PAM)
    except:# PAM on the left
        left_search = SeqUtils.nt_search(seq[:len(PAM)], PAM)
        if len(left_search)>1:
            logger.info('5\' PAM detected')
            PAM_index = left_search[1]
        else:
            right_search = SeqUtils.nt_search(seq[-len(PAM):], PAM)
            if len(right_search)>1:
                PAM_index = len(seq)-len(PAM)
            else:
                logger.info("PAM: %s not found in %s. Set PAM index to 20"%(PAM,seq))
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
        otseq = offtargets[j]['seq']
        bulge_flag = offtargets[j]['bulge_flag']
        dist = offtargets[j]['dist']

        if bulge_flag == False:
            b_dist = levenshtein_two_matrix_rows(otseq, target_seq,PAM)
            k = 0
            line_number += 1
            y = y_offset + line_number * box_size
            for i, (c, r) in enumerate(zip(otseq, target_seq)):
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
        if bulge_flag: #theres a bulge seq present
            k = 0
            line_number += 1
            y = y_offset + line_number * box_size
            for i, (c, r) in enumerate(zip(otseq, realigned_target_seq)):
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

        reads_text = dwg.text(str(seq['reads']), insert=(
            box_size * (len(target_seq) + 1) + 20, y_offset + box_size * (line_number + 2) - 2),
                              fill='black', style="font-size:15px; font-family:Courier")
        dwg.add(reads_text)
        mismatch_text = dwg.text(dist, insert=(
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

        # else:
        #     reads_text = dwg.text(str(seq['reads']), insert=(
        #     box_size * (len(target_seq) + 1) + 20, y_offset + box_size * (line_number + 1) + 5),
        #                           fill='black', style="font-size:15px; font-family:Courier")
        #     dwg.add(reads_text)
        #
        #     mismatch_text = dwg.text(dist, insert=(
        #         box_size * (len(target_seq) + 1) + 130, y_offset + box_size * (line_number + 2) - 2),
        #                              fill='black', style="font-size:15px; font-family:Courier")
        #     dwg.add(mismatch_text)
        #     annot_text = dwg.text(seq['annot'], insert=(
        #         box_size * (len(target_seq) + 1) + 200, y_offset + box_size * (line_number + 1) + 5),
        #                           fill='black', style="font-size:15px; font-family:Courier")
        #     dwg.add(annot_text)
        #     mismatch_text = dwg.text(seq['coord'], insert=(
        #         box_size * (len(target_seq) + 1) + 380, y_offset + box_size * (line_number + 1) + 5),
        #                              fill='black', style="font-size:15px; font-family:Courier")
        #     dwg.add(mismatch_text)
        #
        #     dist_text = dwg.text(dist, insert=(
        #         box_size * (len(target_seq) + 1) + 130, y_offset + box_size * (line_number + 1) - 2),
        #                              fill='black', style="font-size:15px; font-family:Courier")
        #     dwg.add(dist_text)
        #     ## TH added coords
        #
        #     reads_text02 = dwg.text(u"\u007D", insert=(
        #     box_size * (len(target_seq) + 1) + 7, y_offset + box_size * (line_number + 1) + 5),
        #                             fill='black', style="font-size:23px; font-family:Courier")
        #     dwg.add(reads_text02)
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
    parser.add_argument("-tseq", "--target_seq", help="targetseq with PAM", required=True)
    parser.add_argument("--PAM", help="PAM sequence", default="NGG")
    args = parser.parse_args()

    print(args)

    visualizeOfftargets(args.identified_file, args.outfile, args.title, args.target_seq, args.PAM)


if __name__ == "__main__":
    main()

