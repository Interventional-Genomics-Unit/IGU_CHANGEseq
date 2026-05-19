import subprocess
import gzip
import os
import pandas as pd
import yaml
from utility import get_paths,write_yaml_to_file

GeneSource = ['Gnomon', 'Curated Genomic', 'tRNAscan-SE', 'Curated Genomic%2Ccmsearch',
    'BestRefSeq', 'BestRefSeq%2CGnomon', 'cmsearch']
GeneBiotype = [
    'rRNA', 'misc_RNA', 'D_segment', 'ncRNA_pseudogene', 'antisense_RNA', 'ncRNA', 'Y_RNA', 'C_region', 'snRNA',
    'lncRNA', 'miRNA', 'V_segment_pseudogene', 'tRNA', 'scRNA', 'other', 'transcribed_pseudogene', 'V_segment',
    'vault_RNA', 'snoRNA', 'C_region_pseudogene', 'J_segment', 'protein_coding', 'pseudogene', 'RNase_P_RNA',
    'J_segment_pseudogene', 'telomerase_RNA', 'RNase_MRP_RNA']
RNASource = ['Gnomon', 'Curated Genomic', 'tRNAscan-SE', 'BestRefSeq', 'cmsearch']
RNABiotype = ['rRNA', 'antisense_RNA', 'Y_RNA', 'unknown', 'V_gene_segment', 'RNase_P_RNA',
    'snRNA', 'J_gene_segment', 'mRNA', 'miRNA', 'tRNA', 'miRNA_primary_transcript', 'scRNA', 'lnc_RNA_pseudogene',
    'scaRNA', 'C_gene_segment', 'vault_RNA', 'snoRNA', 'RNase_MRP_RNA', 'V_gene_segment_pseudogene', 'pseudogene',
    'lnc_RNA', 'C_gene_segment_pseudogene', 'telomerase_RNA', 'J_gene_segment_pseudogene', 'D_gene_segment',]
RNAExperiment = [
    'COORDINATES: polyA evidence [ECO:0006239]',
    'COORDINATES: cap analysis [ECO:0007248]',
    'COORDINATES: cap analysis [ECO:0007248] and polyA evidence [ECO:0006239]'
]
RNATag =['RefSeq Select', 'RefSeq Plus Clinical']

CDSSource = [
    'Gnomon', 'Curated Genomic', 'BestRefSeq']

ignore_types = {
    "biological_region", "enhancer", "silencer", "transcriptional_cis_regulatory_region",
    "protein_binding_site", "nucleotide_motif", "non_allelic_homologous_recombination_region",
    "recombination_feature", "promoter", "sequence_feature", "meiotic_recombination_region",
    "mobile_genetic_element", "DNaseI_hypersensitive_site", "conserved_region", "origin_of_replication",
    "tandem_repeat", "repeat_instability_region", "mitotic_recombination_region", "enhancer_blocking_element",
    "sequence_alteration", "TATA_box", "region", "response_element", "chromosome_breakpoint",
    "sequence_secondary_structure", "locus_control_region", "matrix_attachment_site",
    "epigenetically_modified_region", "replication_regulatory_region", "direct_repeat", "insulator",
    "minisatellite", "repeat_region", "CAAT_signal", "dispersed_repeat", "microsatellite", "inverted_repeat",
    "nucleotide_cleavage_site", "sequence_comparison", "GC_rich_promoter_region", "replication_start_site",
    "imprinting_control_region", "regulatory_region", "CAGE_cluster", "TSS", "sequence_alteration_artifact",
    "centromere", "match", "cDNA_match", "D_loop"
}
seqids = {"NC_060925.1": "chr1",
        "NC_060926.1": "chr2",
        "NC_060927.1": "chr3",
        "NC_060928.1": "chr4",
        "NC_060929.1": "chr5",
        "NC_060930.1": "chr6",
        "NC_060931.1": "chr7",
        "NC_060932.1": "chr8",
        "NC_060933.1": "chr9",
        "NC_060934.1": "chr10",
        "NC_060935.1": "chr11",
        "NC_060936.1": "chr12",
        "NC_060937.1": "chr13",
        "NC_060938.1": "chr14",
        "NC_060939.1": "chr15",
        "NC_060940.1": "chr16",
        "NC_060941.1": "chr17",
        "NC_060942.1": "chr18",
        "NC_060943.1": "chr19",
        "NC_060944.1": "chr20",
        "NC_060945.1": "chr21",
        "NC_060946.1": "chr22",
        "NC_060947.1": "chrX",
        "NC_060948.1": "chrY",}





##### ------------------------


genes =  {'transcripts': {'exons': {}}
rna = {'exons': {}}

def preprocess_gff(gff_file):
    # gff_file = "/groups/clinical/projects/clinical_shared_data/hg38/annotations/GCF_000001405.40_GRCh38.p14_genomic.gff"
    names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attrs"]

    df = pd.read_csv(
        gff_file, sep='\t', comment="#",
        names=names,
        dtype={
            "seqid": str, "source": str, "type": str, "start": int, "end": int,
            "score": str, "strand": str, "phase": str, "attrs": str
        }).drop(columns=["score", "phase"])

    cur_gene, cur_mrna = None, None
    records = {}
    for i, entry in df.iterrows():
        gff_record = entry.to_dict()

        if entry["type"] == "gene":

        attrs = {}
        for attribute in entry["attrs"].split(";"):
            k, v = attribute.split("=")
            attrs[k] = v
        break


def process_refseq(tmp_output, output):
    '''
    ftp : https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz

    Makes
    A) a bed file from refseqs = 282,614 lines
    B) a bed file of the most recent genes = 66688 lines
    '''

    chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
              '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT', 'Y', 'X']

    #in
    ref_out = gzip.open(output, 'wt')

    labels = ['bin', 'id', 'chrom', 'strand', 'txStart', 'txEnd',
              'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds',
              'score', 'name', 'cdsStartStat', 'cdsEndStat','exonFrames']

    cnt = 0
    for line in gzip.open(tmp_output, 'rt'):
        tokens = line.split('\t')
        tid, chrom, strand,tstart,tend = tokens[1:6]
        cds_start, cds_end = tokens[6], tokens[7]
        exons_start, exon_end = tokens[9], tokens[10]
        gname, frames = tokens[12], tokens[-1].split('\n')[0]

        if tokens[2].replace('chr', "") in chroms:
            tid = tokens[1]
            cnt+=1
            if tid.startswith('X') == False:
                eid = '-'
                line_out = [chrom,tstart,tend,'|'.join([strand,tid,eid,gname,cds_start,cds_end, exons_start,exon_end,frames])]
                ref_out.write('\t'.join(line_out) + '\n')
    ref_out.close()

def get_refseq(ftp_path,tmp_output):
    cmd = "wget " + ftp_path + " -O " +  tmp_output
    subprocess.check_call(cmd, shell=True)

def get_refseq_genome(ftp_path,tmp_output):
    #https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz
    cmd = "wget " + ftp_path + " -O " + tmp_output
    subprocess.check_call(cmd, shell=True)
    cmd = "gunzip " + tmp_output
    subprocess.check_call(cmd, shell=True)



def makefiles(p_dir):
    paths_dict = get_paths(p_dir)

    #file = outdir + "paths.txt" ## for writing paths

    #if reset_output:
    #    print("Reseting changeseqs annotation file to "+ reset_output)
    #    write_path(file, reset_output)
    #else:
    ## get infiles
    ftp_path = paths_dict['ftps']['refseq_txt']
    output = paths_dict['refseq']
    tmp_output = p_dir + "/data/" + str(os.path.basename(ftp_path))

    print("downloading Refseq from " + ftp_path)
    get_refseq(ftp_path, tmp_output)

    print("Cleaning Refseq and converting to bed file")
    process_refseq(tmp_output, output)

    print("Cleaning up tmp files")
    os.remove(tmp_output)


