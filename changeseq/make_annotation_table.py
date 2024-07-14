import subprocess
import gzip
import os

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


def write_path(file,output):
    with open(file, "w") as f:
        f.write(output)
    f.close()


def makefiles(ftp_path,outdir,reset_output,p_dir):
    file = outdir + "\paths.txt" ## for writing paths

    if reset_output:
        print("Reseting changeseqs annotation file to "+ reset_output)
        write_path(file, reset_output)
    else:
        print("downloading Refseq from " + ftp_path)
        tmp_output = outdir + "tmp_ncbiRefSeq.txt.gz"
        output = outdir + "ncbiRefSeq.bed.gz"
        get_refseq(ftp_path, tmp_output)

        print("Cleaning Refseq and converting to bed file")
        process_refseq(tmp_output, output)

        print("Storing file path to ", output)
        write_path(file, output)

        print("Cleaning up tmp files")
        os.remove(tmp_output)


