import os
import pandas as pd
import numpy as np
#matplotlib.use('Agg')  # Use a non-GUI backend
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

global colors
colors = {'sample1': '#8FBC8F', 'sample2': '#029386'}#, 'T': '#D8BFD8', 'C': '#8FBC8F', 'N': '#AFEEEE', 'R': '#3CB371', '-': '#E6E6FA'}

def check_file(file):
    if os.path.isfile(file):
        pass
    else:
        print(file + " is not in the your-CHANGEseq-analysis-path/identified/ directory")

def parse_df(file, threshold = 6):
    df = pd.read_csv(file)
    df = df[df['Nuclease_Read_Count'] >= threshold]

    ## Removing Gaps and merging
    df['Site_Sequence_NoGaps'] = df['Site_Sequence'].fillna(df['Site_Sequence_Gaps_Allowed'].str.replace("-",""))
    df = df[['Genomic Coordinate','Nuclease_Read_Count','Control_Read_Count','Site_Sequence_NoGaps',
        'Site_Sequence','Site_Sequence_Gaps_Allowed','Site_Substitution_Number',
        'DNA_Bulge','RNA_Bulge','Gene_Name','Feature']]

    return df

def join_replicates(sample1_identified_file,sample2_identified_file,threshold = 6):
    sample1_df = parse_df(sample1_identified_file,threshold)
    sample2_df = parse_df(sample2_identified_file,threshold)
    df = sample1_df.merge(sample2_df, on=['Genomic Coordinate','Site_Sequence_NoGaps'], how='outer',suffixes = [".Rep1", ".Rep2"])

    #df = sample1_df.join(sample2_df.set_index('Site_Sequence_NoGaps'),on = 'Site_Sequence_NoGaps',how = 'outer',lsuffix=".Rep1",rsuffix=".Rep2")
    df['Nuclease_Read_Count.Rep1'] = df['Nuclease_Read_Count.Rep1'].fillna(0)
    df['Nuclease_Read_Count.Rep2'] = df['Nuclease_Read_Count.Rep2'].fillna(0)
    df['Gene_Name'] = df['Gene_Name.Rep1'].fillna(df['Gene_Name.Rep2'])
    df['Site_Sequence'] = df['Site_Sequence.Rep1'].fillna(df['Site_Sequence.Rep2'])
    df['Site_Sequence_Gaps_Allowed'] = df['Site_Sequence_Gaps_Allowed.Rep1'].fillna(df['Site_Sequence_Gaps_Allowed.Rep2'])
    #df['Genomic Coordinate'] = df['Genomic Coordinate.Rep1'].fillna(df['Genomic Coordinate.Rep2'])
    df['Site_Substitution_Number'] = df['Site_Substitution_Number.Rep1'].fillna(df['Site_Substitution_Number.Rep2'])
    df['DNA_Bulge'] = df['DNA_Bulge.Rep1'].fillna(df['DNA_Bulge.Rep2'])
    df['RNA_Bulge'] = df['RNA_Bulge.Rep1'].fillna(df['RNA_Bulge.Rep2'])
    df['Feature'] = df['Feature.Rep1'].fillna(df['Feature.Rep2'])

    df['Number of Replicates Sites found'] = [2 if (i * j) > 1 else False for i, j in zip(df['Nuclease_Read_Count.Rep1'], df['Nuclease_Read_Count.Rep2'])]

    df = df[['Genomic Coordinate','Nuclease_Read_Count.Rep1','Nuclease_Read_Count.Rep2',
           'Control_Read_Count.Rep1','Control_Read_Count.Rep2','Site_Sequence_NoGaps',
        'Site_Sequence','Site_Sequence_Gaps_Allowed','Site_Substitution_Number',
        'DNA_Bulge','RNA_Bulge','Gene_Name','Feature']]


    return df


def swarm_plot(df,name,figout):
    sns.set_style('white')
    x1, x2 = list(df['Nuclease_Read_Count.Rep1']), list(df['Nuclease_Read_Count.Rep2'])
    y = [(x + y) / 2 if (x * y) > 0 else x + y for x, y in zip(x1, x2)]
    y = np.log2(np.array(y))
    df['Log2 Mean Reads'] = y
    df['Guide Name'] = [name] * len(df['Log2 Mean Reads'])
    df['Shared'] = [True if (i * j) > 0 else False for i,j in zip(x1,x2)]
    palette=['lightgrey','royalblue']

    # Highlight on target if present

    if len(df.loc[df['Site_Substitution_Number'] + df['RNA_Bulge'].fillna(0) + df['DNA_Bulge'].fillna(0) ==0,'Shared']):
        df.loc[(df['Site_Substitution_Number'] + df['RNA_Bulge'].fillna(0) + df['DNA_Bulge'].fillna(0) == 0), 'Shared'] = 'on target'
        df = df.sort_values('Shared',ascending = False)
        palette = ['#6D091F','royalblue', 'lightgrey']
    plt.figure(figsize = (6,6))
    g = sns.swarmplot(data = df,
                      x= 'Guide Name',
                      y= 'Log2 Mean Reads',
                      hue='Shared',
                      palette=palette
                      )
    plt.yticks(np.arange(np.log2(10), np.log2(100000), step=np.log2(10)), [10, 100, 1000, 10000, 100000])
    plt.title(name)
    plt.ylabel("Read Counts")
    plt.xlabel("")
    #plt.show()
    plt.savefig(figout, bbox_inches='tight')
    plt.close(figout)


def scatter_plot(x1,x2,name,figout):
    x = np.log2(np.array(x1)+1)
    y = np.log2(np.array(x2)+1)

    pearR = np.corrcoef(x1, x2)[1, 0]
    plt.figure(figsize=(4, 4))
    plt.scatter(x[x * y>0], y[x * y>0],color = colors.values()[1],label="rho= %s" % (round(pearR,3)))
    plt.xticks(np.arange(np.log2(10), np.log2(100000), step=np.log2(10)),[10,100,1000,10000,100000])
    plt.yticks(np.arange(np.log2(10), np.log2(100000), step=np.log2(10)),[10,100,1000,10000,100000])
    plt.legend(loc=1)
    plt.title(name)
    #plt.show()
    plt.savefig(figout, bbox_inches='tight')
    plt.close(figout)


def calc_jaccard(Rep1_unique,Rep2_unique,shared):
    total = Rep1_unique + Rep2_unique + shared
    return round(float(shared)/float(total)*100,2)

def create_pseudo_sample(x1,x2):
    pseudo_sample = []
    for i in range(len(x1)):
        pseudo_sample.append((x1[i]+ x2[i])/2)
    return pseudo_sample

def median_normalization(read_counts1,read_counts2):
    # Similar to DESeq2
    #https://divingintogeneticsandgenomics.com/post/details-in-centered-log-ratio-clr-normalization-for-cite-seq-protein-count-data/
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
    joined['Number of Replicates Sites found'] = [2 if (i * j) > 1 else False for i, j in
                                              zip(joined['Nuclease_Read_Count.Rep1'], joined['Nuclease_Read_Count.Rep2'])]
    return joined

def vennplot_replicates(joined,sample1,sample2,figout):
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

    plt.savefig(figout, bbox_inches='tight')
    #plt.show()
    return ja



def repCombiner(sample1,sample2,name,analysis_folder,read_threshold = 6):
    ## Input
    # sample1, sample2  = "781_02471_rep1", "781_02471_rep2"
    # analysis_folder = '/groups/clinical/projects/Assay_Dev/CHANGEseq/CS_12/'
    # name = "NA02471 spCas9 PRF1 sgRNA_781"  # for labeling
    # read_threshold = 6


    ##  Inputs
    sample1_identified_file = analysis_folder + 'identified/' + sample1 + '_identified_matched_annotated.csv'
    sample2_identified_file = analysis_folder + 'identified/' + sample2 + '_identified_matched_annotated.csv'
    check_file(sample1_identified_file)
    check_file(sample2_identified_file)

    ## Outputs
    joined_out = os.path.join(output_folder, sample1) + 'JOINED_RAW' + sample2 + '.csv'
    joined_normalized_out = os.path.join(output_folder,sample1) + 'JOINED_NORMALIZED' + sample2 + '.csv'
    venn_out_without_normalize = analysis_folder + 'visualization/' + name.replace(" ",
                                                                                   "_") + "_replicate_venndiagram.png"
    venn_out = analysis_folder + 'visualization/' + name.replace(" ",
                                                                 "_") + "_replicate_without_normalizing_venndiagram.png"
    scatter_out = analysis_folder + 'visualization/' + name.replace(" ", "_") + "_replicate_scatterplot.png"
    swarm_out = analysis_folder + 'visualization/' + name.replace(" ", "_") + "_replicate_swarmplot.png"

    # Join without normalizing
    print("Joining...")
    print(sample1_identified_file)
    print("with")
    print(sample2_identified_file)

    joined = join_replicates(sample1_identified_file, sample2_identified_file, threshold=read_threshold)
    joined.to_csv(joined_out, index=False)


    print("Writing raw un-normalized vendiagram...")
    print(venn_out_without_normalize)
    sim = vennplot_replicates(joined, sample1, sample2, venn_out_without_normalize)
    print(sim)

    joined_normalized = normalize(joined, threshold=read_threshold)
    joined_normalized.to_csv(joined_normalized_out, index=False)

    # plotting
    sim = vennplot_replicates(joined_normalized, sample1, sample2, venn_out)
    print(sim)

    x1, x2 = list(joined_normalized['Nuclease_Read_Count.Rep1']), list(joined_normalized['Nuclease_Read_Count.Rep2'])

    scatter_plot(x1, x2, name, scatter_out)
    swarm_plot(joined_normalized, name, swarm_out)

'''
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
    
df = pd.read_csv("/groups/clinical/projects/Assay_Dev/CHANGEseq/CASAFE/casoffinder/hg38_CasSAFE_casoffinder.txt",sep = "\t")
df = df.drop(columns = ["Guide_ID","crRNA"])
joined_normalized = joined_normalized.rename(columns = {'Genomic Coordinate':'Coordinates'})

df2 = df.join(joined_normalized.set_index('Coordinates'), on = 'Coordinates')
out = os.path.join(analysis_folder, 'identified',
                                         name) + '_JOINED_NORMALIZED_casoffinder.csv'

df2[~df2['DNA'].isna()]   
    
'''


def parse_args():
    mainParser = argparse.ArgumentParser()
    mainParser.add_argument('--sample1', '-s1', help='name of replicate 1. Must match sample in manifest')
    mainParser.add_argument('--sample2', '-s2', help='name of replicate 2. Must match sample in manifest')
    #mainParser.add_argument('--file1', '-f1', help='absolute path of sample1 identified_annotated.csv file')
    #mainParser.add_argument('--file2', '-f2', help= 'absolute path of sample2 identified_annotated.csv file')
    mainParser.add_argument('--output', '-o', help='output directory')
    mainParser.add_argument('--name', '-n', help='Sample name for labeling graphs')
    mainParser.add_argument('--read_threshold', '-rt', help='limit sample combining to a min read count threshold',default =6)

    return mainParser.parse_args()


def main():
    args = parse_args()
    repCombiner(args.sample1, args.sample2,args.name, args.output, int(args.read_threshold))

if __name__ == '__main__':
    main()