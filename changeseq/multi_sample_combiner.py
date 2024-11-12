import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


global colors
colors = {'sample1': '#8FBC8F',
          'sample2': '#029386'}  # , 'T': '#D8BFD8', 'C': '#8FBC8F', 'N': '#AFEEEE', 'R': '#3CB371', '-': '#E6E6FA'}


def parse_df(file,description,threshold=6):
    df = pd.read_csv(file)
    df = df[df['Nuclease_Read_Count'] >= threshold]

    ## Removing Gaps and merging
    df['Site_Sequence_NoGaps'] = df['Site_Sequence'].fillna(df['Site_Sequence_Gaps_Allowed'].str.replace("-", ""))
    df = df[['Genomic Coordinate', 'Nuclease_Read_Count', 'Control_Read_Count', 'Site_Sequence_NoGaps',
             'Site_Sequence', 'Site_Sequence_Gaps_Allowed', 'Site_Substitution_Number',
             'DNA_Bulge', 'RNA_Bulge', 'Gene_Name', 'Feature']]
    df['Description'] = [str(description)] * len(df['Gene_Name'])

    return df


def join_replicates(sample1_df, sample2_df, suffixes,threshold=6):
    df = sample1_df.merge(sample2_df, on=['Genomic Coordinate', 'Site_Sequence_NoGaps'], how='outer',
                          suffixes=suffixes)

    # df = sample1_df.join(sample2_df.set_index('Site_Sequence_NoGaps'),on = 'Site_Sequence_NoGaps',how = 'outer',lsuffix=".Rep1",rsuffix=".Rep2")
    df['Nuclease_Read_Count'+ suffixes[0]] = df['Nuclease_Read_Count' + suffixes[0]].fillna(0)
    df['Nuclease_Read_Count' + suffixes[1]] = df['Nuclease_Read_Count' + suffixes[1]].fillna(0)
    df['Control_Read_Count' + suffixes[1]] = df['Control_Read_Count' + suffixes[1]].fillna(0)
    df['Control_Read_Count' + suffixes[0]] = df['Control_Read_Count' + suffixes[0]].fillna(0)
    df['Site_Sequence'] = df['Site_Sequence' + suffixes[0]].fillna(df['Site_Sequence' + suffixes[1] ])
    df['Site_Sequence_Gaps_Allowed'] = df['Site_Sequence_Gaps_Allowed' + suffixes[0]].fillna(
        df['Site_Sequence_Gaps_Allowed' + suffixes[1]])
    df['Site_Substitution_Number'] = df['Site_Substitution_Number' + suffixes[0]].fillna(df['Site_Substitution_Number'+suffixes[1]])
    df['DNA_Bulge'] = df['DNA_Bulge' + suffixes[0]].fillna(df['DNA_Bulge'+suffixes[1]])
    df['RNA_Bulge'] = df['RNA_Bulge' + suffixes[0]].fillna(df['RNA_Bulge'+suffixes[1]])
    df['Gene_Name'] = df['Gene_Name' + suffixes[0]].fillna(df['Gene_Name'+suffixes[1]])
    df['Feature'] =  df['Feature' + suffixes[0]].fillna(df['Feature'+suffixes[1]])

    df = df.drop(
        columns=['DNA_Bulge' + suffixes[0],'DNA_Bulge' + suffixes[1],
                 'RNA_Bulge' + suffixes[0], 'RNA_Bulge' + suffixes[1],
                 'Feature' + suffixes[0], 'Feature' + suffixes[1],
                 'Gene_Name' + suffixes[0], 'Gene_Name' + suffixes[1],
                 'Site_Substitution_Number' + suffixes[0],'Site_Substitution_Number' + suffixes[1],
                 'Site_Sequence_Gaps_Allowed' + suffixes[0],'Site_Sequence_Gaps_Allowed' + suffixes[1],
                 'Site_Sequence' + suffixes[0],'Site_Sequence' + suffixes[1]])

    return df


def scatter_by_overlap(x1, x2, name):
    x = np.log2(np.array(x1) + 1)
    y = np.log2(np.array(x2) + 1)
    pearR = np.corrcoef(x1, x2)[1, 0]
    colors = ['#8FBC8F' if j != 0 and k != 0 else 'gray' for j, k in zip(x, y)]
    plt.figure(figsize=(4, 4))

    plt.scatter(x, y, c=colors, label="r= %s" % (round(pearR, 3)))

    # plt.scatter(x[x * y > 0], y[x * y > 0], color=colors.values()[1], label="r= %s" % (round(pearR, 3)))
    plt.xticks(np.arange(np.log2(1), np.log2(100000), step=np.log2(10)), [0, 10, 100, 1000, 10000, 100000])
    plt.yticks(np.arange(np.log2(1), np.log2(100000), step=np.log2(10)), [0, 10, 100, 1000, 10000, 100000])
    plt.legend(loc=1)
    # A = (x + y) / 2
    # M = x - y
    # plt.scatter(x=A, y=M, s=10, c = colors)  # s is point size
    plt.title(name)
    plt.show()


def swarm_plot(df, name, figout):
    x1, x2 = list(df['Nuclease_Read_Count.Rep1']), list(df['Nuclease_Read_Count.Rep2'])
    y = [(x + y) / 2 if (x * y) > 0 else x + y for x, y in zip(x1, x2)]
    y = np.log2(np.array(y))
    df['Log2 Mean Reads'] = y
    df['Guide Name'] = [name] * len(df['Log2 Mean Reads'])
    df['Shared'] = [True if (i * j) > 0 else False for i, j in zip(x1, x2)]
    plt.figure(figsize=(5, 5))
    g = sns.swarmplot(data=df,
                      x='Guide Name',
                      y='Log2 Mean Reads',
                      hue='Shared',
                      palette=['lightgrey', 'royalblue'])
    plt.yticks(np.arange(np.log2(10), np.log2(100000), step=np.log2(10)), [10, 100, 1000, 10000, 100000])
    plt.title(name)
    plt.ylabel("Read Counts")
    plt.xlabel("")
    plt.show()
    plt.savefig(figout, bbox_inches='tight')
    plt.close(figout)



def calc_jaccard(Rep1_unique, Rep2_unique, shared):
    total = Rep1_unique + Rep2_unique + shared
    return round(float(shared) / float(total) * 100, 2)

def geo_mean(iterable):
    a = np.array(iterable)
    return a.prod()**(1.0/len(a))


def create_pseudo_sample(counts):
    pseudo_sample = []
    for i in range(counts.shape[0]):
        pseudo_sample.append(geo_mean(counts[i,:]))
    return np.array(pseudo_sample)


def median_normalization(counts):
    # https://divingintogeneticsandgenomics.com/post/details-in-centered-log-ratio-clr-normalization-for-cite-seq-protein-count-data/
    pseudo_sample = create_pseudo_sample(counts) # (log_cnt1,log_cnt2)

    norm_count = pseudo_sample[:, None] / counts
    norm_count[np.isnan(norm_count)] = 0
    size_factor = np.zeros(np.shape(norm_count[0,:]))
    for i in range(norm_count.shape[1]):
        x = norm_count[:,i]
        x = x[~np.isnan(x)]
        size_factor[i] = np.median(x[x>0])
    new_counts = counts * size_factor
    #new_read_counts1 = [round(x / size_factor1, 0) for x in read_counts1]
    #new_read_counts2 = [round(x / size_factor2, 0) for x in read_counts2]
    new_counts = new_counts.round(0)
    return new_counts

def normalize(joined):
    counts = np.array(joined.loc[:,joined.columns.str.startswith('Nuc')].fillna(0))
    joined.loc[:, joined.columns.str.startswith('Nuc')] = median_normalization(counts)
    return joined

def multi_combine(file,name,analysis_folder, read_threshold=6):
    # analysis folde must have all '_identified_matched_annotated.csv' files present
    # name = "spCas9 CASAFE"
    # analysis_folder = '/groups/clinical/projects/Assay_Dev/CHANGEseq/CASAFE/'
    # file =  analysis_folder + "multiple_combine.csv"
    # read_threshold = 6

    joined_out = analysis_folder + name.replace(' ',"") +'_multi_sample_combine.csv'
    joined_norm = analysis_folder + name.replace(' ', "") + '_multi_sample_normalized_combine.csv'

    first_file = False
    with open(file,"r") as f:
        for line in f:
            line = line.strip().split(",")
            if line[0].startswith('sample'):
                first_file = True
                pass
            else:
                id_file = analysis_folder + 'identified/' + line[0] + '_identified_matched_annotated.csv'
                df = parse_df(id_file, line[1])
                df = df.rename(columns={'Nuclease_Read_Count': 'Nuclease_Read_Count.' + line[0],
                                                 'Control_Read_Count':'Control_Read_Count.' +line[0],
                                        'Description':'Description.' +line[0]})
                if first_file:
                    joined = df
                    previous_suffix = '.' + line[0]
                    first_file = False
                else:
                    suffix = '.' + line[0]
                    joined = join_replicates(joined, df, suffixes = [previous_suffix,suffix],threshold = read_threshold)
                    previous_suffix = suffix

    joined_normalized = normalize(joined)
    joined.to_csv(joined_out,index = False)
    joined_normalized = joined_normalized[[u'Genomic Coordinate', u'Nuclease_Read_Count.spCas9_12878_CASAFE_r1',
     u'Nuclease_Read_Count.spCas9_12878_CASAFE_r2',
    u'Nuclease_Read_Count.spCas9_24631_CASAFE_r1',
     u'Nuclease_Read_Count.spCas9_24631_CASAFE_r2',
    u'Site_Sequence_NoGaps',u'Site_Substitution_Number',
    u'DNA_Bulge', u'RNA_Bulge',
    u'Gene_Name', u'Feature',
    u'Site_Sequence_Gaps_Allowed']]
    joined_normalized.to_csv(joined_norm,index = False)



def parse_args():
    mainParser = argparse.ArgumentParser()
    mainParser.add_argument('--file', '-f', help='manifest file of samples to be combined')
    mainParser.add_argument('--output', '-o', help='output directory')
    mainParser.add_argument('--name', '-n', help='Sample name for labeling graphs')
    mainParser.add_argument('--read_threshold', '-rt', help='limit sample combining to a min read count threshold',
                            default=6)

    return mainParser.parse_args()


def main():
    args = parse_args()
    try:
        multi_combine(args.file, args.name, args.output, args.read_threshold)
    except Exception as e:
        print('Error combined replicates')