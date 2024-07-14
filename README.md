# IGU CHANGE-seq
 An extension of CHANGE-seq for Off-Target Site Nomination
 https://github.com/tsailabSJ/changeseq
 
IGU Added output

* manifest.yaml file creation by .csv input
* Optional UMI deduplication
* replicate combination
* normalization of replicates
* Fastq QC using Fastp
* Identified site bam file and alignment historgram
* annotation to sites
* annotations to visualization 


 ## Installation
 git clone the IGU_CHANGEseq into a designated folder using the following command
 
 ```
 git clone https://github.com/Interventional-Genomics-Unit/IGU_CHANGEseq
 ```
 
 Create a conda environment and import all required packages in one line below:
 
 ```
 conda env create -f cseq.yaml
 conda activate cseq
 ```

In order to set up the annotation table for annotating OT sites run makefiles to download RefSeq annotations. Additional commands can change the location and FTP link (--outdir --ftp_path)

```
python /home/thudson/projects/IGU_CHANGEseq/changeseq/changeseq.py makefiles

```
 
# Usage

The change-seq pipeline requires a manifest yaml file specifying input files, output directory, and pipeline parameters. this file is created using two csv files. One is specfic to paramters and the other are the input paramters. Once the csv files are created, users can simply run the change-seq wrapper 

```
cd IGU_CHANGEseq/changeseq/
bash CHANGEseq_wrapper.sh test/parameters.csv test/manifest.csv all all

```


## Example Output


# Writing A Manifest File
When running the end-to-end analysis functionality of the CHANGEseq package a number of inputs are required. To simplify the formatting of these inputs and to encourage reproducibility, these parameters are inputted into the pipeline via a manifest. This first manifest.csv sample takes into the 


    - For each sample, you must provide the following parameters:
        - `target`: Target sequence for that sample. Accepts degenerate bases.
        - `read1`: The absolute path to the .FASTQ(.gz) file containing the read1 reads.
        - `read2`: The absolute path to the .FASTQ(.gz) file containing the read2 reads.
        - `controlread1`: The absolute path to the .FASTQ(.gz) file containing the control read1 reads.
        - `controlread2`: The absolute path to the .FASTQ(.gz) file containing the control read2 reads.
        - `description`: A brief description of the sample
     

 # Writing A Parameters File
 This file is a .csv with the first column the paramter name and the second column the paraemeter name. See IGU_CHANGEseq/changeseq/test/parameters.csv

- `reference_genome`: The absolute path to the reference genome FASTA file.
- `analysis_folder`: The absolute path to the folder in which all pipeline outputs will be saved.
- 'raw_fastq_folder':The absolute path to the folder in which all raw fastq files are stored.
- `bwa`: The absolute path to the `bwa` executable
- `samtools`: The absolute path to the `samtools` executable
- `read_threshold`: The minimum number of reads at a location for that location to be called as a site. We recommend leaving it to the default value of 4.
- `window_size`: Size of the sliding window, we recommend leaving it to the default value of 3.
- `mapq_threshold`: Minimum read mapping quality score. We recommend leaving it to the default value of 50.
- `start_threshold`: Tolerance for breakpoint location. We recommend leaving it to the default value of 1.
- `gap_threshold`: Distance between breakpoints. We recommend leaving it to the default value of 3 for Cas9.
- `mismatch_threshold`: Number of tolerated gaps in the fuzzy target search setp. We recommend leaving it to the default value of 6.
- `read_length`: Fastq file read length, default is 151.
- `PAM`: PAM sequence, default is NGG.
- `genome`: used for homer peak annotation, e.g., hg19, hg38, mm9, or mm10.
- `merged_analysis`: Whether or not the paired read merging step should takingTrue
- 'dedup_umi': Whether or not the dedupliucation of UMIs step should takingTrue or False
- 'bc_pattern': If umi deduplication is taking place, which barcode pattern should be given to UMI tools

# Pipeline Output
When running the full pipeline, the results of each step are outputted to the `output_folder` in a separate folder for each step. The output folders and their respective contents are as follows:

- `output_folder/aligned`: Contains an alignment `.sam`, alignment `.bam`, sorted `bam`, and `.bai` index file for each sample.
- `output_folder/fastq`: Merged `.fastq.gz` files for each sample.
- `output_folder/identified`: Contains tab-delimited `.txt` files for each sample containing the identified DSBs, control DSBs, filtered DSBs, and read quantification.
- `output_folder/visualization`: Contains a `.svg` vector image representing an alignment of all detected off-targets to the targetsite for each sample.
