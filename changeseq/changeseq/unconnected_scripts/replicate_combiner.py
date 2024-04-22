import os
import pandas as pd
import math
import numpy as np
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

global colors
colors = {'sample1': '#8FBC8F', 'sample2': '#029386'}#, 'T': '#D8BFD8', 'C': '#8FBC8F', 'N': '#AFEEEE', 'R': '#3CB371', '-': '#E6E6FA'}

def parse_df(file, threshold = 6):
    df = pd.read_csv(file)
    df = df[df['Nuclease_Read_Count'] >= threshold]
    df = df[['Genomic Coordinate','Nuclease_Read_Count','Control_Read_Count',
        'Site_Sequence','Site_Sequence_Gaps_Allowed','Site_Substitution_Number',
        'DNA_Bulge','RNA_Bulge','Gene_Name','Feature']]
    return df

def join_replicates(sample1_identified_file,sample2_identified_file,threshold = 6):
    sample1_df = parse_df(sample1_identified_file,threshold)
    sample2_df = parse_df(sample2_identified_file,threshold)
    df = sample1_df.join(sample2_df.set_index('Genomic Coordinate'),on = 'Genomic Coordinate',how = 'outer',lsuffix=".Rep1",rsuffix=".Rep2")
    df['Nuclease_Read_Count.Rep1'] = df['Nuclease_Read_Count.Rep1'].fillna(0)
    df['Nuclease_Read_Count.Rep2'] = df['Nuclease_Read_Count.Rep2'].fillna(0)

    print('Number of Rep1 sites:', len(df[df['Nuclease_Read_Count.Rep2']==0]))
    print('Number of Rep2 sites:', len(df[df['Nuclease_Read_Count.Rep1'] == 0]))
    print('Number of shared Sites:', len(df[df['Nuclease_Read_Count.Rep2'] * df['Nuclease_Read_Count.Rep1'] > 0]))
    print('Number of different Sites:', len(df[df['Nuclease_Read_Count.Rep2'] * df['Nuclease_Read_Count.Rep1'] == 0]))
    print('Number in both', len(df['Nuclease_Read_Count.Rep1']))

    return df

def scatter_plot(x1,x2,name):
    x = np.log2(np.array(x1)+1)
    y = np.log2(np.array(x2)+1)

    pearR = np.corrcoef(x1, x2)[1, 0]
    plt.scatter(x[x * y>0], y[x * y>0],color = colors.values()[1],label="r= %s" % (round(pearR,3)))
    plt.xticks(np.arange(np.log2(10), np.log2(100000), step=np.log2(10)),[10,100,1000,10000,100000])
    plt.yticks(np.arange(np.log2(10), np.log2(100000), step=np.log2(10)),[10,100,1000,10000,100000])
    plt.legend(loc=1)
    plt.title(name)
    plt.show()

def create_pseudo_sample(x1,x2):
    pseudo_sample = []
    for i in range(len(x1)):
        pseudo_sample.append((x1[i]+ x2[i])/2)
    return pseudo_sample

def calc_jaccard(Rep1_unique,Rep2_unique,shared):
    total = Rep1_unique + Rep2_unique + shared
    return round(float(shared)/float(total)*100,2)



def median_normalization(read_counts1,read_counts2):
    #https://divingintogeneticsandgenomics.com/post/details-in-centered-log-ratio-clr-normalization-for-cite-seq-protein-count-data/
    log_cnt1 = [np.log(c+1) for c in read_counts1]
    log_cnt2 = [np.log(c+1) for c in read_counts2]
    pseudo_sample = create_pseudo_sample(read_counts1,read_counts2)#(log_cnt1,log_cnt2)
    norm_count1 = [y / x for x,y in zip(pseudo_sample,read_counts1)]
    norm_count2 = [y / x for x,y in zip(pseudo_sample,read_counts2)]
    size_factor1 = np.median(norm_count1[norm_count1>0])
    size_factor2 = np.median(norm_count2[norm_count2>0])
    new_read_counts1 = [round(x / size_factor1,0)for x in read_counts1]
    new_read_counts2 = [round(x / size_factor2,0) for x in read_counts2]
    #scatter_plot(new_read_counts1, new_read_counts2)
    return new_read_counts1, new_read_counts2

def normalize(joined,threshold=6):
    read_counts1, read_counts2 = list(joined['Nuclease_Read_Count.Rep1']), list(joined['Nuclease_Read_Count.Rep2'])
    rep1,rep2 = median_normalization(read_counts1, read_counts2)
    joined['Nuclease_Read_Count.Rep1'] = [x if x >=threshold else 0 for x in rep1]
    joined['Nuclease_Read_Count.Rep2'] = [x if x >= threshold else 0 for x in rep2]
    joined = joined[(joined['Nuclease_Read_Count.Rep1'] + joined['Nuclease_Read_Count.Rep2'] !=0)]
    return joined

def vennplot_replicates(joined,sample1,sample2):
    Rep1_unique = len(joined[joined['Nuclease_Read_Count.Rep2'] == 0])
    Rep2_unique =len(joined[joined['Nuclease_Read_Count.Rep1'] == 0])
    shared = len(joined[joined['Nuclease_Read_Count.Rep2'] * joined['Nuclease_Read_Count.Rep1'] > 0])
    ja= calc_jaccard(Rep1_unique,Rep2_unique,shared)
    values = (Rep1_unique, Rep2_unique,shared)
    names = (sample1,sample2)
    plt.figure(figsize=(4, 4))
    v =venn2(subsets = values,set_labels=names,set_colors=(colors['sample1'],colors['sample2']),alpha = 0.5)
    v.get_label_by_id("A").set_fontsize(8)
    v.get_label_by_id("B").set_fontsize(8)
    v.get_label_by_id("A").set_y(0.6)
    v.get_label_by_id("B").set_y(0.6)
    v.get_label_by_id("A").set_x(len(sample1)/100.0-0.4)
    v.get_label_by_id("B").set_x(len(sample1)/100.0)
    plt.annotate("% Replicate Sites Overlap " + str(ja), xy=v.get_label_by_id('010').get_position() +
                                           np.array([0, -0.5]), xytext=(-60, -30), ha='center',
                 textcoords='offset points')
                 #bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
                 #arrowprops=dict(arrowstyle='->',
                 #                connectionstyle='arc', color='gray'))
    plt.show()
    return ja




# x = ['/groups/clinical/projects/Assay_Dev/IGU_CHANGEseq/CS_03/identified/spCas9_G6PC1_10584_rep1_identified_matched.txt', '//groups/clinical/projects/Assay_Dev/IGU_CHANGEseq/CS_03/identified/spCas9_G6PC1_10584_rep2_identified_matched.txt']
sample1 = "spCas9_CASAFE2_dual_guide_rep2"
sample2 = "spCas9_CASAFE2_dg_cs5_r2"
sample1_name,sample2_name = "spCas9 CaSAFE2 Rep1", "spCas9 CaSAFE2 Rep2"
name = "spCas9 CASAFE2 Between Runs"
r = 6
analysis_folder = ('/groups/clinical/projects/Assay_Dev/IGU_CHANGEseq/CASAFE/')
sample1_identified_file = os.path.join(analysis_folder, 'identified', sample1) + '_identified_matched_annotated.csv'
analysis_folder = ('/groups/clinical/projects/Assay_Dev/IGU_CHANGEseq/CS_05/12878/')
sample2_identified_file = os.path.join(analysis_folder, 'identified', sample2) + '_identified_matched_annotated.csv'

joined = join_replicates(sample1_identified_file,sample2_identified_file,threshold = r)
joined_out = os.path.join(analysis_folder, 'identified', sample1) + 'JOINED' + sample2 + '.csv'
sim =vennplot_replicates(joined, sample1_name, sample2_name)
print(sim)
x1,x2 = list(joined['Nuclease_Read_Count.Rep1']),list(joined['Nuclease_Read_Count.Rep2'])
scatter_plot(x1,x2,name = name)

joined.to_csv(joined_out,index = False)


data_similarity = []
read_thresholds = [6,12,24,48,96,192]
normalize_data = True
for r in read_thresholds:
    joined = join_replicates(sample1_identified_file,sample2_identified_file,threshold = r)
    if normalize_data:
        joined_normalized = normalize(joined,threshold=r)
        sim =vennplot_replicates(joined_normalized, sample1_name, sample2_name)
    else:
        sim = vennplot_replicates(joined,sample1,sample2)
    data_similarity.append(sim)


#joined_out = os.path.join(analysis_folder, 'identified', sample1) + 'JOINED' + sample2 + '.csv'
#joined.to_csv(joined_out,index = False)

