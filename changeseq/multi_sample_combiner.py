import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import logging
from matplotlib_venn import venn2
from visualize import *
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

#############
## author: Taylor Hudson
## Innovative Genomics Institute, UC Berkeley
##
## script combines change-seq repilcates and normalizes read counts based on
## user selected method (median-ratio, rpm). Alignment plots and tables are remade on normalized counts
## along with strip plots, venn diagrams and swarmlots.
################


#global colors
#colors = {'sample1': '#8FBC8F',
#          'sample2': '#029386'}  # , 'T': '#D8BFD8', 'C': '#8FBC8F', 'N': '#AFEEEE', 'R': '#3CB371', '-': '#E6E6FA'}
logger = logging.getLogger('root')
logger.propagate = False

global colors
colors = ['#406061', '#96BA66', '#B45C20','#718CA0',"#5C842F"]

def check_file(file):
    if os.path.isfile(file):
        pass
    else:
        logger.info(f"{file} is not your directory")

def get_read_depth(qcfiles):
    depths = []
    for file in qcfiles:
        check_file(file)
        line = open(file).readline()
        if 'total_reads:' in line:
            depths.append(int(line.split(':')[-1].strip().replace(",","")))
        else:
            logger.info(f'Total read depth for {file} could not be found for tpm normalization')
            logger.info('Re-run QC or choose another normalization method')
            exit()
    return depths

def find_mm_and_bulge(df):
    subs,insertions,deletions = [],[],[]
    for i,data in df.iterrows():
        db = 0
        rb = 0
        mm=0
        try:
            mm = int(data['Site_Substitution_Number'])
        except ValueError:
            for i,j in zip(data['Site_Sequence_Gaps_Allowed'][:-3],data['Aligned_Target_Sequence'][:-3]):
                if i==j:
                    pass
                elif i == "-":
                    rb += 1
                elif j =='-':
                    db += 1
                else:
                    mm+=1
        subs.append(mm)
        insertions.append(rb)
        deletions.append(db)
    return subs,insertions,deletions

def choose_cols(df):

    keep_cols = ['Genomic Coordinate']
    keep_cols += [x for x in df.columns if 'Nuclease_Read_Count' in x or 'Control_Read_Count' in x]

    keep_cols += ['Site_Sequence','Aligned_Site_Sequence','Site_Substitution_Number',
                  'Target_Sequence', 'Aligned_Target_Sequence',
                  'DNA_Bulge', 'RNA_Bulge','Distance','Number of Replicates Sites found','Percent Total Reads',
                  'Gene_Name', 'Feature', 'Cell', 'MappingPositionStart',
                  'MappingPositionEnd', 'WindowSequence']

    if 'Nuclease Base Edited' in df.columns:
        keep_cols +=['Nuclease Base Edited', 'Nuclease Noise Base Edited','Nuclease pval','Control pval']
    if 'gnomAD.constraint' in df.columns:
        keep_cols += ['gnomAD.constraint', 'HPA.disease_involvement', 'COSMICS.cancer_role',
                      'COSMICS.tier']

    return keep_cols



def parse_df(df,sample):
    ## Filling in Aligned and ungapped sequence
    df.loc[:,'Aligned_Site_Sequence']= df['Site_Sequence_Gaps_Allowed'].fillna(df['Site_Sequence'])
    df.loc[:,'Site_Sequence'] = df['Site_Sequence'].fillna(df['Site_Sequence_Gaps_Allowed'].astype('str').str.replace("-", ""))
    df.loc[:, 'Aligned_Target_Sequence'] = df['Realigned_Target_Sequence'].astype('str').replace('none', np.nan)
    df.loc[:,'Aligned_Target_Sequence'] = df.loc[:,'Aligned_Target_Sequence'].fillna(df['Target_Sequence'])
    df.loc[:,'RNA_Bulge'] = df['RNA_Bulge'].fillna(0)
    df.loc[:,'DNA_Bulge'] = df['DNA_Bulge'].fillna(0)
    df.loc[:,'Genomic Coordinate'] = df['Genomic Coordinate'].str.split('-',expand=True)[0] + df['Strand']

    dist = df.loc[:,['Site_Substitution_Number', 'RNA_Bulge', 'DNA_Bulge']].sum(1)
    i = list(df.columns).index('DNA_Bulge')
    df.insert(i+1,'Distance',dist)

    df.insert(i+2,'Number of Replicates Sites found',1)

    x = df.loc[:, 'Nuclease_Read_Count'] / df.loc[:, 'Nuclease_Read_Count'].sum()
    df.insert(i + 3,'Percent Total Reads',x.round(5) *100)
    keep_cols = choose_cols(df)
    df = df.loc[:,keep_cols]


    df = df.rename(columns={'Nuclease_Read_Count': 'Nuclease_Read_Count.' + sample,
                            'Control_Read_Count': 'Control_Read_Count.' + sample,
                            'Description':'Description.' +sample})

    return df


def join_replicates(sample1_df, sample2_df, suffixes):
    df = sample1_df.merge(sample2_df, on=['Genomic Coordinate', 'Site_Sequence'], how='outer',
                          suffixes=suffixes)


    keep_cols = choose_cols(df)

    for col in df.columns:
        if "Read_Count" in col:
            df.loc[:,col] = df[col].fillna(0)

    for col in keep_cols[keep_cols.index( 'Aligned_Site_Sequence'):]:
        df.loc[:,col] = df[f'{col}{suffixes[0]}'].fillna(df[f'{col}{suffixes[1]}'])
        df = df.drop(columns=[f'{col}{suffixes[0]}', f'{col}{suffixes[1]}'])

    nuclease_read_count_cols = [col for col in df.columns if "Nuclease_Read_Count" in col]
    df.loc[:,'Number of Replicates Sites found'] = (df[nuclease_read_count_cols] > 0).sum(1)


    x = df.loc[:, nuclease_read_count_cols].sum(1) / df.loc[:, nuclease_read_count_cols].sum(1).sum()
    df.loc[:,'Percent Total Reads'] = x.round(5) *100


    return df



def swarm_plot(joined_normalized, name, figout):
    sns.set_style('white')
    columns = [col for col in joined_normalized.columns if col.startswith('Nuclease_Read_Count')]
    count_dict = joined_normalized.loc[:,joined_normalized.columns.str.startswith('Nuclease_Read_Count')].to_dict('list')

    # Finding the mean reading excluding 0 sites
    mean = []
    counts = np.array(list(count_dict.values()))
    for i in range(len(counts[0,:])):
        count = counts[:,i]
        x = sum(count[count>0])
        mean.append(x/sum(count>0))


    df = joined_normalized.copy()
    df.loc[:,'Log2 Mean Reads'] = np.log2(np.array(mean))
    df.loc[:,'Guide Name'] = [name] * len(df['Log2 Mean Reads'])
    df.loc[:,'Shared'] = [f"found in {x} samples" for x in df['Number of Replicates Sites found']]

    # Highlight on target if present
    ontarget_row = df[['Site_Substitution_Number', 'RNA_Bulge', 'DNA_Bulge']].sum(1)==0
    df.loc[ontarget_row , 'Shared'] = 'on target'
    df = df.sort_values('Shared', ascending=False)
    colors2 = ['#96C04B', '#7ACBC8', '#7c98d3', '#97FFFF']
    palette = colors2[0:len(columns)+1]
    plt.figure(figsize=(6, 6))
    g = sns.swarmplot(data=df,
                      x='Guide Name',
                      y='Log2 Mean Reads',
                      hue='Shared',
                      palette=palette
                      )
    plt.yticks(np.arange(np.log2(10), np.log2(100000), step=np.log2(10)), [10, 100, 1000, 10000, 100000])
    plt.title(name)
    plt.ylabel("log2 read counts")
    plt.xlabel("")
    legend = g.get_legend()
    if legend:
        legend.set_title("")
    # plt.show()
    plt.savefig(figout)
    plt.close(figout)

def calc_jaccard(Rep1_unique,Rep2_unique,shared):
    total = Rep1_unique + Rep2_unique + shared
    return round(float(shared)/float(total)*100,2)


def geo_mean(iterable):
    a = np.array(iterable)
    return a.prod()**(1.0/len(a))


def create_pseudo_sample(counts):
    pseudo_sample = []
    for i in range(counts.shape[0]):
        pseudo_sample.append(geo_mean(counts[i,:]))
    return np.array(pseudo_sample)

def median_normalization(count_dict):
    # https://divingintogeneticsandgenomics.com/post/details-in-centered-log-ratio-clr-normalization-for-cite-seq-protein-count-data/
    normalized_data={}

    counts = np.array(list(count_dict.values())).T.astype(float)
    pseudo_sample = create_pseudo_sample(counts) # (log_cnt1,log_cnt2)
    mask = pseudo_sample > 0
    counts_f = counts[mask, :]
    pseudo_f = pseudo_sample[mask]
    norm_count = pseudo_f[:, None] / counts_f
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
    c=0
    for k in count_dict.keys():
        normalized_data[k] = list(new_counts[:,c].astype('int'))
        c+=1

    return normalized_data

def rpm(count_dict,depths):
    # uses a total count scaling method

    min_depth = np.min(depths)

    scaling_factors = {
        sample: min_depth / depths[i]
        for i,sample in enumerate(list(count_dict.keys()))
    }

    normalized_data = {
        sample: (np.array(counts) * scaling_factors[sample]).round(0)
        for sample, counts in count_dict.items()
    }
    return normalized_data


def scatter_plot(x1,x2,name,figout):
    x = np.log2(np.array(x1)+1)
    y = np.log2(np.array(x2)+1)

    pearR = np.corrcoef(x1, x2)[1, 0]
    plt.figure(figsize=(4, 4))
    plt.scatter(x[x * y>0], y[x * y>0],color = colors[3],label="rho= %s" % (round(pearR,3)))
    plt.xticks(np.arange(np.log2(10), np.log2(100000), step=np.log2(10)),[10,100,1000,10000,100000])
    plt.yticks(np.arange(np.log2(10), np.log2(100000), step=np.log2(10)),[10,100,1000,10000,100000])
    plt.legend(loc=1)
    plt.title(name)
    #plt.show()
    plt.savefig(figout, bbox_inches='tight')
    plt.close(figout)

def vennplot_replicates(x1,x2,sample1,sample2,figout):
    font_size = 8
    x1, x2 =np.array(x1),np.array(x2)
    Rep1_unique,Rep2_unique  = sum(x1>0),sum(x2 > 0)
    shared = sum(x1 * x2 > 0)
    ja= calc_jaccard(Rep1_unique,Rep2_unique,shared)
    values = (Rep1_unique, Rep2_unique,shared)

    plt.figure(figsize=(4, 4))
    v =venn2(subsets = values,
             set_labels=(sample1,sample2),
             set_colors=(colors[1],colors[4]),alpha = 0.5)
    for l in ["A","B"]:
        v.get_label_by_id(l).set_fontsize(font_size)
        v.get_label_by_id(l).set_y(0.6)
        v.get_label_by_id(l).set_x(len(sample1)/100.0-0.4)

    plt.annotate("% replicate rites overlap " + str(ja), xy=v.get_label_by_id('010').get_position() +
                                           np.array([0, -0.5]), xytext=(-60, -30), ha='center',
                 textcoords='offset points')
    plt.savefig(figout, bbox_inches='tight')
    #plt.show()
    return ja
def clean_normalized(joined_normalized,read_threshold):
    read_columns = joined_normalized.columns.str.startswith('Nuclease_Read_Count')

    # remove sites below threshold post normalization
    joined_normalized.loc[:, read_columns] = joined_normalized.loc[:,read_columns].applymap(
        lambda x: 0 if x and x < read_threshold else x)
    keep_rows = joined_normalized.loc[:, read_columns].sum(1) != 0
    joined_normalized = joined_normalized.loc[keep_rows, :].copy()
    joined_normalized.loc[:,'Number of Replicates Sites found'] = (joined_normalized.loc[:, read_columns] > 0).sum(1)
    joined_normalized.loc[:,'Percent Total Reads'] =joined_normalized.loc[:, read_columns].sum(1) / joined_normalized.loc[:, read_columns].sum(1).sum()
    joined_normalized.loc[:, 'Percent Total Reads'] =  joined_normalized['Percent Total Reads'].round(5) *100

    # make a simplified version
    keep_cols = choose_cols(joined_normalized)

    joined_simplified_report =  joined_normalized.loc[:, keep_cols ].copy()
    return joined_normalized, joined_simplified_report

# def drop_large_controls(df,read_threshold):
#     df2 = df.copy()
#     df2['Nuclease_Read_Count'] = df2['Nuclease_Read_Count'] - (df2['Control_Read_Count'] * 2)
#     df2 = df2.loc[df2['Nuclease_Read_Count']>=read_threshold,:]
#     df2['Nuclease_Read_Count'] = df2['Nuclease_Read_Count'].round(0).astype('int')
#
#     return df2



def normalize(joined,qcfiles,normalization_method,read_threshold):
    count_dict = joined.loc[:,joined.columns.str.startswith('Nuclease_Read_Count')].to_dict('list')
    ## check if at least 5% of the sites are matching. if not then skip normalization:
    pt = 1-sum([x.count(0) for x in count_dict.values()])/len(list(count_dict.values())[0])

    if pt < 0.05:
        norm_count_dict = count_dict
        logger.info(f'Not enough overlapping sites to normalize {count_dict.keys()}')
    elif normalization_method == "median":
        norm_count_dict = median_normalization(count_dict)
    elif normalization_method == "rpm":
        depths =get_read_depth(qcfiles)
        norm_count_dict =rpm(count_dict,depths)
    else:
        norm_count_dict = count_dict
    for k,v in norm_count_dict.items():
        joined.loc[:, k] = v
    joined_normalized, simplified_report = clean_normalized(joined,read_threshold)

    return joined_normalized,simplified_report

def make_offtarget_dict(joined_normalized,subset):
    sample_df = joined_normalized.loc[joined_normalized[subset]>0].reset_index().copy()
    offtargets = []
    total_seq = sample_df.shape[0]
    target_seq = joined_normalized['Target_Sequence'].iloc[0].replace("-", "")

    for i,row in sample_df.iterrows():
        offtarget_reads = row[subset]
        annot = ""
        if int(row['RNA_Bulge']) + int(row['DNA_Bulge']) == 0:
            no_bulge_offtarget_sequence = row['Site_Sequence']
            bulge_offtarget_sequence = ""
            realigned_target_seq = row['Aligned_Target_Sequence']

        else:
            no_bulge_offtarget_sequence = ""
            realigned_target_seq =  row['Aligned_Target_Sequence']
            bulge_offtarget_sequence = row['Aligned_Site_Sequence']

        coord = row['Genomic Coordinate']

        try:
            if "intergenic" not in row['Feature']:
                annot = row['Gene_Name'] + "," + row['Feature'].replace("non-coding RNA","ncRNA")  # gene name and feature

        except:
            pass

        if no_bulge_offtarget_sequence != '' or bulge_offtarget_sequence != '':
            if no_bulge_offtarget_sequence:
                total_seq += 1
            if bulge_offtarget_sequence:
                total_seq += 1
            offtargets.append({'seq': no_bulge_offtarget_sequence.strip(),
                               'bulged_seq': bulge_offtarget_sequence.strip(),
                               'reads': int(offtarget_reads),
                               'target_seq': target_seq.strip(),
                               'coord': coord,
                               'annot': str(annot),
                               'realigned_target_seq': realigned_target_seq.strip()
                               })
    offtargets = sorted(offtargets, key=lambda x: x['reads'], reverse=True)
    return offtargets, target_seq, total_seq

def process_results(rep_group_name,replicates,infiles,qcfiles,outfolder, normalization_method,read_threshold,PAM):
    '''
    joined and normalizes replicates as indicated in manifest
    '''
    logger.info(f'Joining and normalizing {rep_group_name} replicates')
    processed_outfile =outfolder + "/tables/"+ rep_group_name +'_joined.csv'
    simplified_report_outfile = processed_outfile.replace('.csv','_simplified.csv')
    swarm_plot_out = outfolder + "/visualization/"+ rep_group_name + "_postprocess_swarmplot.png"

    first_file = True
    for i,infile in enumerate(infiles):
        sample = replicates['sample_name'][i]
        df = pd.read_csv(infile)
        df = parse_df(df,sample)
        #df = drop_large_controls(df,read_threshold)

        if first_file:
            joined = df
            previous_suffix = '.' + sample
            first_file = False
        else:
            suffix = '.' +sample
            joined = join_replicates(joined, df, suffixes = [previous_suffix,suffix])
            previous_suffix = suffix


    joined_normalized, simplified_report = normalize(joined,qcfiles,normalization_method,read_threshold)
    joined_normalized.to_csv(processed_outfile, index = False)
    simplified_report.to_csv(simplified_report_outfile, index = False)

    ## plotting
    swarm_plot(joined_normalized, rep_group_name, swarm_plot_out)

    for sample in replicates['sample_name']:
        offtargets, target_seq, total_seq = make_offtarget_dict(joined_normalized,subset='Nuclease_Read_Count.' + sample)
        alignment_plot = outfolder +"/visualization/"+ sample.replace(" ", "_") + "_postprocess_alignment_plot.svg"
        draw_plot(target_seq, offtargets, total_seq, outfile=alignment_plot, title=sample, PAM=PAM)

    for i in range(len(replicates['sample_name'])-1):
        sample_1= replicates['sample_name'][i]  #'Nuclease_Read_Count.' + sample
        for j in range(1,len(replicates['sample_name'])):
            sample_2 = replicates['sample_name'][j]
            x1, x2 = list(joined_normalized[f'Nuclease_Read_Count.{sample_1}']), list(
                joined_normalized[f'Nuclease_Read_Count.{sample_2}'])
            scatter_out = f"{outfolder}/visualization/{sample_1}_&_{sample_2}_postprocess_scatterplot.png"
            venn_out = f"{outfolder}/visualization/{sample_1}_&_{sample_2}_postprocess_venn.png"

            scatter_plot(x1, x2,f"{sample_1} & {sample_2}", scatter_out)
            sim = vennplot_replicates(x1,x2,sample_1, sample_2, venn_out)


def parse_args():
    mainParser = argparse.ArgumentParser()
    mainParser.add_argument('--replicate_sample_names', '-r', help='replicate sample names as indicated on manifest, comma seperated ex. CCR5_rep1,CCR5_rep2')
    mainParser.add_argument('--output', '-o', help='output directory')
    mainParser.add_argument('--name', '-n', help='Sample name for labeling graphs')
    mainParser.add_argument('--infiles', '-f', help='list of i_identified_matched_annotated.csv files in raw_output folder ')
    mainParser.add_argument('--normalization_method', '-nm',
                            help='')
    mainParser.add_argument('--analysis_folder', '-p',
                            help='')
    mainParser.add_argument('--PAM', '-p',
                            help='')
    mainParser.add_argument('--read_threshold', '-rt', help='limit sample combining to a min read count threshold',
                            default=6)

    return mainParser.parse_args()


def main():
    args = parse_args()
    try:
        replicates = {'sample_name': args.replicate_sample_names.split(",")}
        process_results(args.name,replicates,args.infiles,args.qcfiles,args.analysis_folder, args.normalization_method,args.read_threshold,args.PAM)
        process_results(args.file, args.name, args.output, args.read_threshold)
    except Exception as e:
        print('Error combined replicates')
