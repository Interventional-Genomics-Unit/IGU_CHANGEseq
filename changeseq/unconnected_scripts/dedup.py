from collections import defaultdict
import gzip


def deduplicate_fastq(input_file, output_file):
    sequences = defaultdict(list)

    with gzip.open(input_file, 'rt') as f_in, gzip.open(output_file, 'w') as f_out:
        while True:
            header = f_in.readline().strip()
            if not header:  # Check for end of file
                break
            seq = f_in.readline().strip()
            plus = f_in.readline().strip()  # Optional line, typically starts with '+'
            qual = f_in.readline().strip()

            # Store sequence in a dictionary
            sequences[seq].append((header, qual))

        for seq, headers_quals in sequences.items():
            for header, qual in headers_quals:
                f_out.write(header + '\n')
                f_out.write(seq + '\n')
                f_out.write('+\n')
                f_out.write(qual + '\n')


# Usage example

analysis_folder = '/groups/clinical/projects/Assay_Dev/IGU_CHANGEseq/CASAFE/'
sample = 'control_spCas9_CASAFE_dual_guide_rep1'
input_file = "{0}fastq/{1}_merged.fastq.gz".format(analysis_folder,sample)
output_file = "{0}fastq/{1}_dedup_merged.fastq.gz".format(analysis_folder,sample)

deduplicate_fastq(input_file, output_file)