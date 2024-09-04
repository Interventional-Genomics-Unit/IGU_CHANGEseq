from __future__ import print_function
import subprocess
import sys
import os
import argparse
import regex
import re
import HTSeq
import pyfaidx
from findCleavageSites import get_sequence, regexFromSequence, alignSequences, reverseComplement, extendedPattern, realignedSequences


"""
Run samtools:mpileup and get all identified variants in the window sequences
"""

def snpCall(matched_file, reference, bam_file, out, search_radius):
    basename = os.path.basename(out)
    output_folder = os.path.dirname(out)

    # open matched file
    regions = list()
    with open(matched_file, 'rU') as f:
        f.readline()
        for line in f:
            site = line.strip().split('\t')
            #  chromosome, windowStart, windowEnd, strand, bam, region_basename (=Targetsite_Name)
            regions.append([site[0], int(site[-4]) - search_radius, int(site[-3]) + search_radius, '*', bam_file, site[15]])#'_'.join([site[26], site[3]])])

    print('Running samtools:mpileup for %s' % basename) #, file=sys.stderr)
    out_bcf = os.path.join(output_folder, basename + '_output_bcftools')
    if os.path.exists(out_bcf):
        subprocess.check_call('rm -r %s' % out_bcf, shell=True, env=os.environ.copy())
    os.makedirs(out_bcf)
    process_mpileup = open(os.path.join(out_bcf, 'logFile_mpileup'), 'w')

    """
    print('Running samtools:mpileup for %s' % basename)  # , file=sys.stderr)
    out_bcf = os.path.join(output_folder, basename + '_output_bcftools')
    process_mpileup = open(os.path.join(out_bcf, 'logFile_mpileup'), 'w')
    output = os.path.join(out_bcf,  + basename + '_BCFcall.vcf.gz')
    cl_vcf = "bcftools mpileup --fasta-ref %s %s | bcftools call -cv -Oz -o %s" % (reference, bam_file, output)
    subprocess.check_call(cl_vcf, shell=True, env=os.environ.copy(), stderr=process_mpileup, stdout=process_mpileup)
    process_mpileup.close()


    print('Collecting significant variant calls for %s' % basename) #, file=sys.stderr)
    out_svc = os.path.join(output_folder, basename + '_output_svc')
    process_svc = open(os.path.join(out_svc, 'logFile_svc'), 'w')
    output = os.path.join(out_svc, name + '_SIGNFcall.txt')
    cl_sed = "bcftools filter -i 'QUAL>20' %s | sed -n '/##/!p' | awk 'FNR>1' > %s" % (os.path.join(out_bcf, arch), output)
    subprocess.check_call(cl_sed, shell=True, env=os.environ.copy(), stderr=process_svc, stdout=process_svc)
    process_svc.close()

    """
    for item in regions:
        chromosome, windowStart, windowEnd, strand, bam_file, region_basename = item
        region = '%s%s%s%s%s' % (chromosome, ":", int(windowStart), "-", int(windowEnd))
        output = os.path.join(out_bcf, region_basename + '_BCFcall.vcf.gz')

        cl_vcf = "bcftools mpileup --region %s --fasta-ref %s %s | bcftools call -cv -Oz -o %s" % (region, reference, bam_file, output)
        subprocess.check_call(cl_vcf, shell=True, env=os.environ.copy(), stderr=process_mpileup, stdout=process_mpileup)
    process_mpileup.close()

    print('Collecting significant variant calls for %s' % basename) #, file=sys.stderr)
    out_svc = os.path.join(output_folder, basename + '_output_svc')
    if os.path.exists(out_svc):
        subprocess.check_call('rm -r %s' % out_svc, shell=True, env=os.environ.copy())
    os.makedirs(out_svc)
    process_svc = open(os.path.join(out_svc, 'logFile_svc'), 'w')

    bcf_files = [f for f in os.listdir(out_bcf) if os.path.isfile(os.path.join(out_bcf, f))]
    for arch in bcf_files:
        if not arch.startswith('.') and arch.endswith('.vcf.gz'):
            name = arch[:-12]
            output = os.path.join(out_svc, name + '_SIGNFcall.txt')
            cl_sed = "bcftools filter -i 'QUAL>20' %s | sed -n '/##/!p' | awk 'FNR>1' > %s" % (os.path.join(out_bcf, arch), output)
            #cl_sed = "zcat %s | sed -n '/##/!p' | awk 'FNR>1' > %s" % (os.path.join(out_bcf, arch), output)
            subprocess.check_call(cl_sed, shell=True, env=os.environ.copy(), stderr=process_svc, stdout=process_svc)
    process_svc.close()

    print('Consolidating all the significant variant calls for %s' % basename) #, file=sys.stderr)
    header = ['targetsite', 'site_name', 'chromosome', 'one_based_position', 'reference', 'variant', 'quality', 'genotype', 'depth', 'PL']
    variants = list()

    svc_files = [f for f in os.listdir(out_svc) if os.path.isfile(os.path.join(out_svc, f))]
    for arch in svc_files:
        if not arch.startswith('.') and arch.endswith('.txt'):

            # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
            ##ALT=<ID=*,Description="Represents allele(s) other than observed.">
            ##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
            ##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of raw reads supporting an indel">
            ##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of raw reads supporting an indel">
            ##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
            ##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
            ##INFO=<ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
            ##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
            ##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
            ##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">
            ##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
            ##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
            ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            ##INFO=<ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)">
            ##INFO=<ID=AF2,Number=1,Type=Float,Description="Max-likelihood estimate of the first and second group ALT allele frequency (assuming HWE)">
            ##INFO=<ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count (no HWE assumption)">
            ##INFO=<ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">
            ##INFO=<ID=FQ,Number=1,Type=Float,Description="Phred probability of all samples being the same">
            ##INFO=<ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">
            ##INFO=<ID=G3,Number=3,Type=Float,Description="ML estimate of genotype frequencies">
            ##INFO=<ID=HWE,Number=1,Type=Float,Description="Chi^2 based HWE test P-value based on G3">
            ##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
            ##bcftools_callVersion=1.11+htslib-1.11

            tag = arch[:-14]
            f = open(os.path.join(out_svc, arch), 'r')
            reads = f.readlines()
            f.close()

            for line in reads:
                item = line.split()
                if 'INDEL' in item[7]:
                    variants.append(
                        [basename, tag] + item[:2] + item[3:6] + [str(int(item[9][0])) + '|' + str(int(item[9][2]))] +
                        [item[7].split(';')[3][3:]]+ ['_'.join(item[9][4:].split(','))])
                else:
                    variants.append(
                        [basename, tag] + item[:2] + item[3:6] + [str(int(item[9][0])) + '|' + str(int(item[9][2]))] +
                        [item[7].split(';')[0][3:]]  + ['_'.join(item[9][4:].split(','))])

    out_file = open(out + '_mpileupCall.txt', 'w')
    print('\t'.join(header), file=out_file)
    for item in variants:
        print('\t'.join(item),file=out_file)
    out_file.close()

    print('Cleaning up directive for %s' % basename, file=sys.stderr)
    #subprocess.check_call('rm -r %s' % out_vcf, shell=True, env=os.environ.copy())
    subprocess.check_call('rm -r %s' % out_bcf, shell=True, env=os.environ.copy())
    subprocess.check_call('rm -r %s' % out_svc, shell=True, env=os.environ.copy())

    print('Done running samtools:mpileup for %s' % basename) #, file=sys.stderr)
    return variants


"""
Obtain variant off-target sequences
"""
def realignVariantBulge(bulge_sequence, window_sequence_variant, bulge_strand):
    bseq = bulge_sequence.replace('-', '')
    if bulge_strand == '+':
        m_bulge = re.search(bseq, window_sequence_variant, re.I)
    else:
        m_bulge = re.search(bseq, reverseComplement(window_sequence_variant), re.I)
    variant_bseq = m_bulge.group()
    variant_bseq = variant_bseq[:bulge_sequence.find('-')] + '-' + variant_bseq[bulge_sequence.find('-'):]
    return variant_bseq


def SNPreader(snp_file):
    ga = HTSeq.GenomicArray("auto", stranded=False, typecode='O')

    for snp in snp_file:
        basename, snpID, chromosome, one_based_position, reference, variant, quality, genotype, depth, PL = snp
        position = int(one_based_position) - 1
        key = '_'.join([basename, chromosome])
        ga[HTSeq.GenomicInterval(chromosome, position, position + 1, ".")] = [position, reference, variant, genotype, quality, key, depth]
    return ga


def arrayOffTargets(matched_file, search_radius):
    offtargets_dict = {}
    gi_dict = {}

    with open(matched_file, 'r') as g:
        g.readline()
        for line in g:
            site = line.strip().split('\t')

            Chromosome = site[0]
            start = int(site[-4]) - search_radius
            end = int(site[-3]) + search_radius
            Name = site[15]

            offtargets_dict[Name] = site

            #  create a genomic interval for each window sequence
            gi_dict[Name] = HTSeq.GenomicInterval(Chromosome, start, end, ".")
    return offtargets_dict, gi_dict


def snpAdjustment(matched_file, snp_file, out, mismatch_threshold, search_radius):
    output_file = open(out + '_Variants.txt', 'w')
    print('#Chromosome', 'Start', 'End', 'Genomic Coordinate', 'Nuclease_Read_Count', 'Strand',
          'Variant_WindowSequence',
          'Variant_Site_SubstitutionsOnly.Sequence', 'Variant_Site_SubstitutionsOnly.NumSubstitutions',
          'Variant_Site_SubstitutionsOnly.Strand',
          'Variant_Site_GapsAllowed.Sequence', 'Variant_Site_GapsAllowed.Length', 
          'Variant_Site_GapsAllowed.Substitutions', 'Variant_Site_GapsAllowed.Insertions', 'Variant_Site_GapsAllowed.Deletions',
          'Variant_Site_GapsAllowed.Strand',
          'Cell', 'Targetsite', 'TargetSequence', 'Variant_RealignedTargetSequence',
          'Reference', 'Variant', 'Genotype', 'Quality', 'Depth',
          sep='\t', file=output_file)
    output_file.close()

    basename = os.path.basename(out)
    offtargets, gi_offtargets = arrayOffTargets(matched_file, search_radius)
    ga_snp = SNPreader(snp_file)

    for name in offtargets:
        variant_flag = False
        site = offtargets[name]
        gi = gi_offtargets[name]

        chromosome = site[0]
        window_sequence = site[-1]
        window_sequence = window_sequence.upper()
        cell, targetsite = site[13], site[14]
        TargetSequence = site[16]
        output01 = site[0:6]
        output03 = [cell, targetsite, TargetSequence]
        ots_nb, ots_bu = site[7], site[9]

        #  obtain variant window sequence
        wkey = '_'.join([basename, chromosome])
        insert_start, insert_end, insert_var, snp_data = list(), list(), list(), {}

        for i, v in ga_snp[gi].steps():
            if v:

                position, reference, variant, genotype, quality, key, depth = v
                if key == wkey:
                    variant = variant.split(',')[0]
                    for n, pos in enumerate(range(gi.start, gi.end)):
                        if pos == int(position):
                            insert_var.append(variant.lower())
                            insert_start.append(n)
                            end_pos = n + len(reference)
                            insert_end.append(end_pos)
                            snp_data[str(position)] = [position, reference, variant, genotype, quality, depth]

        tri = 0
        window_sequence_variant = ''
        for i in range(len(insert_var)):
            variant = insert_var[i]
            pos = insert_start[i]
            window_sequence_variant += window_sequence[tri:pos] + variant.lower()
            tri = insert_end[i]
        window_sequence_variant += window_sequence[tri:]

        #  variant off-target sequences: only proceed if there is a variant in the window sequence
        window_sequence_var = window_sequence_variant.upper()
        if window_sequence_var != window_sequence:
            offtarget_sequence_no_bulge, mismatches, offtarget_sequence_length, chosen_alignment_strand_m, start_no_bulge, end_no_bulge, \
            realigned_target, \
            bulged_offtarget_sequence, length, score, substitutions, insertions, deletions, chosen_alignment_strand_b, bulged_start, bulged_end = \
                alignSequences(TargetSequence, window_sequence_var, max_score=mismatch_threshold)

            variant_ots_no_bulge, variant_ots_bulge = '', ''

            #  get variant sequence if the off-target sequences have changed by considering the variant window
            if ots_nb != offtarget_sequence_no_bulge:
                variant_flag = True
                if chosen_alignment_strand_m == '+':
                    m_no_bulge = re.search(offtarget_sequence_no_bulge, window_sequence_variant, re.I)
                else:
                    m_no_bulge = re.search(offtarget_sequence_no_bulge, reverseComplement(window_sequence_variant), re.I)
                variant_ots_no_bulge = m_no_bulge.group()

            if ots_bu != bulged_offtarget_sequence:
                variant_flag = True
                variant_ots_bulge = realignVariantBulge(bulged_offtarget_sequence, window_sequence_variant, chosen_alignment_strand_b)

            # collect and write variant data if we have variant off-target sequence(s)
            if variant_flag:
                total_genotype, total_reference, total_variant, total_quality = '', '', '', ''
                for pos in snp_data:
                    position, reference, variant, genotype, quality, depth = snp_data[pos]
                    if total_genotype != '':
                        total_genotype += ''.join([':', genotype])
                        total_reference += ''.join([':', reference])
                        total_variant += ''.join([':', variant])
                        total_quality += ''.join([':', quality])
                    else:
                        total_genotype += ''.join([genotype])
                        total_reference += ''.join([reference])
                        total_variant += ''.join([variant])
                        total_quality += ''.join([quality])

                output02 = [variant_ots_no_bulge, mismatches, chosen_alignment_strand_m,
                            variant_ots_bulge, length, substitutions, insertions, deletions, chosen_alignment_strand_b]
                output04 = [total_reference, total_variant, total_genotype, total_quality, depth]
                output_line = output01 + [window_sequence_variant] + output02 + output03 + [realigned_target] + output04

                with open(out + '_Variants.txt', 'a') as output_file:
                    print(*output_line, sep='\t', file=output_file)


"""
Main function
"""
def getVariants(matched_file, ref, bam_file, out, search_radius, mismatch_threshold):
    basename = os.path.basename(out)
    output_folder = os.path.dirname(out)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    snp_file = snpCall(matched_file, ref, bam_file, out, search_radius)

    print('Obtaining Variant Off-Target Sequences for %s' % basename, file=sys.stderr)
    snpAdjustment(matched_file, snp_file, out, mismatch_threshold, search_radius)


def main():
    parser = argparse.ArgumentParser(description='Implement samtools:mpileup to identify genomic variants and adjust the off-target sequence when required.')
    parser.add_argument('--matched_file', help="full_path_to/matched file in 'identified' folder", required=True)
    parser.add_argument('--ref', help="Reference Genome Fasta", required=True)
    parser.add_argument('--bam', help="Sorted BAM file", required=True)
    parser.add_argument('--search_radius', help="Search radius around the position window", default=20, type=int)
    parser.add_argument('--mismatch_threshold', help='Maximum score threshold', default=7, type=int)
    parser.add_argument('--out', help="Output file basename, with full path", required=True)
    args = parser.parse_args()

    getVariants(args.matched_file, args.ref, args.bam, args.out, args.search_radius, args.mismatch_threshold)

if __name__ == "__main__":
    main()


# df = pd.read_csv("/groups/clinical/projects/Assay_Dev/CHANGEseq/CS_09/unmerged/hg38/variants/spCas9_471_ADA_1545_R2_Variants.txt",sep= "\t")