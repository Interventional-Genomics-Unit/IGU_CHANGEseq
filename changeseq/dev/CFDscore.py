import pickle
import re




def load_model_params(score_name, models_dir):
    try:
        if score_name == 'cfd':
            mm_scores = pickle.load(open(models_dir + '/mismatch_score.pkl', 'rb'))
            pam_scores = pickle.load(open(models_dir + '/PAM_scores.pkl', 'rb'))
            return mm_scores, pam_scores

    except Exception:
        print("File containing scores could not be opened")


def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)


def cfd_score(seq1, seq2, models_dir):
    # Doench 2016 off-target scoring
    # Doench, Fusi, et al.  Nature Biotechnology 34, 184–191 2016

    mm_scores, pam_scores = load_model_params('cfd', models_dir)
    pam = seq2[-3:]
    seq1 = seq1.upper().replace('T', 'U')
    seq2 = seq2[:-3].upper().replace('T', 'U')
    m_seq1 = re.search('[^ATCGU\\-]', seq1)
    m_seq2 = re.search('[^ATCGU\\-]', seq2)

    score =1

    if (m_seq1 is None) and (m_seq2 is None):
        if seq1 != seq2:
            shorter, longer = sorted([seq1, seq2], key=len)
            for i in range(-len(shorter), 0):
                if (seq1[i] != seq2[i]):
                    key = 'r' + seq1[i] + ':d' + revcom(seq2[i]) + ',' + str(20 + i + 1)
                    score *= mm_scores[key]

            score *= pam_scores[pam[-2:]]
    else:
        score = -1
    return round(score,4)



def cfd_spec_score(sum_cfd_scores):

    #on_target_seq is spacer site and off_target includes pam


    guide_cfd_score = 100 / (100+sum_cfd_scores)
    guide_cfd_score = round(guide_cfd_score*100,2)
    return guide_cfd_score


def score(input,output,score_name,models_dir):
    #input = '/home/thudson/projects/CrisprScoringHub/test/azimuth_input.txt'
    score_dict = {'name': [],
                  'grna_seq':[],
                  'context_seq':[]}

    with open(input,"r") as inp:
        for line in inp:
            if line.startswith('name')==False:
                name,grna_seq,context_seq = line.strip().split("\t")
                score_dict['name'].append(name)
                score_dict['grna_seq'].append(grna_seq)
                score_dict['context_seq'].append(context_seq)

    if score_name == 'cfd':
        scores = []
        for seqs in zip(score_dict['grna_seq'], score_dict['context_seq']):
            seq1,seq2 = seqs
            scores.append(cfd_score(seq1[:-3], seq2, models_dir))

    score_dict['score'] = scores

    with open(output,'w') as out:
        out.write('\t'.join(score_dict.keys()) + '\n')
        for i in range(len(score_dict['name'])):
            line =[]
            for v in score_dict.values():
                line.append(str(v[i]))
            out.write('\t'.join(line) +'\n')

