from subprocess import Popen, PIPE
import os
from os import listdir, remove
from collections import defaultdict
from annotate import Transcript
from CFDscore import cfd_score

'''
Created by T. Hudson using scripts\software from;

Cas-OFFinder 3.0.0 beta (Jan 24 2021)
Copyright (c) 2021 Jeongbin Park and Sangsu Bae
Website: http://github.com/snugel/cas-offinder

Usage: cas-offinder {input_filename|-} {C|G|A}[device_id(s)] {output_filename|-}
(C: using CPUs, G: using GPUs, A: using accelerators)
'''


# INPUT LEG
def make_casoffinder_input(infile, fasta_fname, pam, pamISfirst, guidelen, grna, gid, casoff_params):
    ## create input file for cas-offinder
    mm, RNAbb, DNAbb, PU = casoff_params

    with open(infile, 'w') as f:
        f.writelines(fasta_fname + "\n")
        line = 'N' * guidelen
        line = line + pam + " "+str(DNAbb) +" "+str(RNAbb) + "\n"
        f.writelines(line)

        dpam = 'N' * len(pam)
        line = grna + dpam + " " + str(mm) + " " + gid + "\n"
        f.writelines(line)


# INPUT/OUTPUT LEG
def cas_offinder_bulge(input_filename, output_filename, cas_off_expath, bulge):
    '''
     The cas-offinder off-line package contains a bug that doesn't allow bulges
    This script is partially a wrapper for cas-offinder to fix this bug
     created by...
    https://github.com/hyugel/cas-offinder-bulge
    '''
    # INPUT LEG
    fnhead = input_filename.replace("_input.txt", "")
    id_dict = {}
    if bulge == True:
        with open(input_filename) as f:
            chrom_path = f.readline()
            pattern, bulge_dna, bulge_rna = f.readline().strip().split()
            isreversed = False  ## PAM 5
            bulge_dna, bulge_rna = int(bulge_dna), int(bulge_rna)
            targets = [line.strip().split() for line in f]
            rnabulge_dic = defaultdict(lambda: [])
            bg_tgts = defaultdict(lambda: set())
            for raw_target, mismatch, gid in targets:
                target = raw_target.rstrip('N')
                len_pam = len(raw_target) - len(target)
                bg_tgts['N' * bulge_dna + target + 'N' * len_pam].add(mismatch)
                id_dict['N' * bulge_dna + target + 'N' * len_pam] = gid
                for bulge_size in range(1, bulge_dna+1):
                    for i in range(1, len(target)):
                        bg_tgt = 'N' * (bulge_dna - bulge_size) + target[:i] + 'N' * bulge_size + target[i:] + 'N' * len_pam
                        bg_tgts[bg_tgt].add(mismatch)
                        id_dict[bg_tgt] = gid

                for bulge_size in range(1, bulge_rna+1):
                    for i in range(1, len(target)-bulge_size):
                        bg_tgt = 'N' * (bulge_dna + bulge_size) + target[:i] + target[i+bulge_size:] + 'N' * len_pam
                        bg_tgts[bg_tgt].add(mismatch)
                        rnabulge_dic[bg_tgt].append((i, int(mismatch), target[i:i+bulge_size]))
                        id_dict[bg_tgt] = gid
            seq_pam = pattern[-len_pam:]


        with open(fnhead + '_bulge.txt', 'w') as f:
            f.write(chrom_path)
            f.write(bulge_dna*'N' + pattern + '\n')
            cnt = 0
            for tgt, mismatch in bg_tgts.items():
                f.write(tgt + ' ' + str(max(mismatch)) + ' ' + '\n')
                cnt+=1

        # THIS FILE PATH IS SUPPLIED TO CASOFF-FINDER
        casin = fnhead + '_bulge.txt'

    else:
        nobulge_dict = {}
        with open(input_filename) as inf:
            for line in inf:
                entry = line.strip().split(' ')
                if len(entry) > 2 and len(entry[-1]) > 3:
                    seq, mm, gid = entry
                    nobulge_dict[seq] = [gid, mm]
        casin = input_filename

    print("Created temporary file (%s)." % (casin))
    # THIS FILE PATH IS SUPPLIED TO CASOFF-FINDER
    outfn = fnhead + '_temp.txt'
    print("Running Cas-OFFinder (output file: %s)..." % outfn)

    p = Popen([cas_off_expath, casin, 'C', outfn])
    ret = p.wait()
    if ret != 0:
        print("Cas-OFFinder process was interrupted!")
        exit(ret)
    print("Processing output file...")

    # OUTPUT LEG
    with open(outfn) as fi, open(output_filename, 'w') as fo:
        fo.write('Coordinates\tDirection\tGuide_ID\tBulge type\tcrRNA\tDNA\tMismatches\tBulge Size\n')\
        #fo.write('Guide_ID\tBulge type\tcrRNA\tDNA\tChromosome\tPosition\tDirection\tMismatches\tBulge Size\n')
        ot_coords = []
        for line in fi:
            entries = line.strip().split('\t')
            ncnt = 0

            if bulge == False:
                gid, mm = nobulge_dict[entries[0]]
                coord = entries[1] + ':' +entries[2] + '-' +str(int(entries[2]) + len(entries[0]))
                to_write = '\t'.join([coord,entries[4],gid,'X',entries[0],entries[3],str(entries[5]),'0'])

                fo.write(to_write+'\n')
                ot_coords.append(coord)
            else:
                for c in entries[0]:
                    if c == 'N':
                        ncnt += 1
                    else:
                        break

                if entries[0] in rnabulge_dic:
                    gid = id_dict[entries[0]]
                    for pos, query_mismatch, seq in rnabulge_dic[entries[0]]:
                        tgt = (entries[0][ncnt:ncnt + pos] + seq + entries[0][ncnt + pos:-len_pam] + seq_pam,
                                   entries[3][ncnt:ncnt + pos] + '-' * len(seq) + entries[3][ncnt + pos:])


                        if query_mismatch >= int(entries[5]):

                            start = int(entries[2]) + (ncnt if (not isreversed and entries[4] == "+") or (isreversed and ncnt > 0 and entries[4] == "-") else 0)
                            coord = entries[1] + ':' + str(start) + '-' + str(int(start) + len(tgt[1]))
                            to_write = '\t'.join([coord, entries[4], gid, 'RNA', tgt[0], tgt[1], str(entries[5]), str(len(seq))])

                            ot_coords.append(coord)
                            fo.write(to_write+'\n')

                else:
                    gid = id_dict[entries[0]]
                    nbulge = 0
                    if isreversed:
                        for c in entries[0][:-ncnt][len_pam:]:
                            if c == 'N':
                                nbulge += 1
                            elif nbulge != 0:
                                break
                        tgt = (seq_pam + entries[0][:-ncnt][len_pam:].replace('N', '-'), entries[3][:-ncnt])
                    else:
                        for c in entries[0][ncnt:][:-len_pam]:
                            if c == 'N':
                                nbulge += 1
                            elif nbulge != 0:
                                break
                        tgt = (entries[0][ncnt:][:-len_pam].replace('N', '-') + seq_pam, entries[3][ncnt:])
                    start =int(entries[2]) + (ncnt if (not isreversed and entries[4] == "+") or (isreversed and ncnt > 0 and entries[4] == "-") else 0)
                    btype = 'X' if nbulge == 0 else 'DNA'
                    coord = entries[1] + ':' + str(start) + '-' + str(int(start) + len(tgt[1]))
                    ot_coords.append(coord)
                    to_write = '\t'.join([coord, entries[4], gid, btype, tgt[0], tgt[1], str(entries[5]), str(nbulge)])
                    fo.write(to_write +'\n')
        remove(fnhead + '_temp.txt')
        remove(casin)



def score_ot(crrna,otseq,models_dir):
    score = 0
    if '-' not in crrna[:-3]:
        if '-' not in otseq[:-4]:

            score = cfd_score(crrna[:-3].upper(), otseq.upper(),models_dir)
    return score


def annotate_ots(output_filename,annote_path,models_dir):
    '''
    Reads output, Scores each off-target seq and annotates each off_target
    '''

    editor = 'spCas9'
    coords = []
    scores = []
    spec_scores = {}
    res = open(output_filename, 'r').readlines()
    for line in res[1:]:
        line = line.strip().split('\t')
        coords.append(line[0])
        if line[2] not in spec_scores.keys():
            spec_scores[line[2]] = 0 if editor == 'spCas9' else '.'
        if editor == 'spCas9':
            score = score_ot(line[4], line[5],models_dir)
            if score > 0 and score != 1:
                spec_scores[line[2]] = spec_scores[line[2]] + score
        else:
            score = '.'
        scores.append(score)

    Transcript.load_transcripts(annote_path,coords)

    new_lines = []
    annotate_out = output_filename.replace('_output', '_annotated')
    with open(annotate_out, 'w') as anout:
        anout.write(res[0].strip() + '\tAnnotation\tScore\n')
        cnt = 0
        for line in res[1:]:
            line = line.strip().split('\t')
            ans = Transcript.transcript(line[0])

            if ans != 'intergenic':
                tid, gname = ans.tx_info()[1:3]
                feature = ans.feature
                if ans.feature.startswith('interg'):
                    x = '|'.join([feature,tid])
                else:
                    x = '|'.join([feature, gname, tid])
                new_line = line + [x, str(scores[cnt])]
                new_lines.append(new_line)

            else:
                x = 'intergenic'
                new_line = line + [x, str(scores[cnt])]
                new_lines.append(new_line)
            cnt += 1
            anout.write('\t'.join(new_line) + '\n')

    # remove(output_filename)

def run_casoffinder(resultsfolder,
                    gseq,
                    fasta_fname,
                    gid,
                    cas_off_expath,
                    genome_name,
                    casoff_params,
                    annote_path,
                    models_dir):

    guidelen = len(gseq)
    editor = 'spCas9'
    pamISfirst = False
    pam = 'NGG'

    if casoff_params[1:3] == (0, 0):
        bulge = False
    else:
        bulge = True
    # for each editor type find off_targets

    # make input file
    infile = resultsfolder + genome_name + "_" + editor + "_" + gid + "_casoffinder_input.txt"

    make_casoffinder_input(infile, fasta_fname, pam, pamISfirst, guidelen, gseq.upper(), gid, casoff_params)

    output_filename = infile.replace('_input.txt', '_output.txt')

    cas_offinder_bulge(infile, output_filename, cas_off_expath, bulge)

    annotate_ots(output_filename,annote_path,models_dir)



p_dir = os.path.dirname(os.path.realpath(__file__))
resultsfolder = '/groups/clinical/projects/Assay_Dev/CHANGEseq/CASAFE/'
#models_dir = '/home/thudson/projects/editability/src/pkl/models/'
models_dir = p_dir + "/changeseq/dev/"
cas_off_expath = '/home/thudson/miniconda3/envs/edit/bin/cas-offinder'
fasta_fname = '/groups/clinical/projects/clinical_shared_data/hg38/hg38.fa'
annote_path = open(p_dir + "/changeseq/data/paths.txt", "r").readlines()[0]
genome_name = 'hg38'
guide_src_name = 'hg38'
mm, RNAbb, DNAbb, PU = 4, 1, 1, "C"
casoff_params = (mm, RNAbb, DNAbb, PU)
#gseq = "GGAACCCAGCGAGTGAAAGANGG"
gseq = "GGGAACCCAGCGAGTGAAGA"
gid = "CASAFE"

run_casoffinder(resultsfolder,
                    gseq,
                    fasta_fname,
                    gid,
                    cas_off_expath,
                    genome_name,
                    casoff_params,
                    annote_path,
                    models_dir)