# IGU CHANGE-seq
 An extension of CHANGE-seq for Off-Target Site Nomination
 https://github.com/tsailabSJ/changeseq
 
IGU Added Features:

* updated to python3
* manifest.yaml file creation by .csv input
* removal of ME adapters and illumina adapter (unmerged version only)
* Optional UMI deduplication
* replicate combination
* normalization of replicates
* Fastq QC using Fastp
* Identified site bam file and alignment historgram
* annotation to sites
* annotations to visualization
* Fixed variable analysis bug


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

The Tsai lab change-seq pipeline requires a manifest yaml file specifying input files, output directory, and pipeline parameters. In the IGU version this file is created using two csv files 1) A Parameters file 2) A list of Samples and Controls file. Once the csv files are created, users can simply run the change-seq wrapper as shown below. 


The bash wrapper is initiated with the parameters.csv, manifest.csv, command, and sample_name

```
cd IGU_CHANGEseq/changeseq
bash path-to/CHANGEseq_wrapper.sh test/merged_parameters.csv test/manifest.csv all all
```

# Writing A Manifest File
When running the end-to-end analysis functionality of the CHANGEseq package a number of inputs are required. To simplify the formatting of these inputs and to encourage reproducibility, these parameters are inputted into the pipeline via a manifest. This first manifest.csv sample takes into the 


- `sequencing_sample_name`: sample name as indicated in the raw fastq (Leave out the ("_S[#]_R1_001")
- `sample_name`: The sample name as you wish to see on the output file names  
- `control_sequencing_sample_name`: No RNP control sample name as indicated in the raw fastq (Leave out the ("_S[#]_R1_001")
- `target`: Target sequence for that sample. Accepts degenerate bases.
- `description`: A brief description of the sample. Typically the cell or gDNA Source name
     

 # Writing A Parameters File
 This file is a .csv with the first column the paramter name and the second column the paraemeter name. See IGU_CHANGEseq/changeseq/test/parameters.csv

- `reference_genome`: The absolute path to the reference genome FASTA file.
- `analysis_folder`: The absolute path to the folder in which all pipeline outputs will be saved.
- `raw_fastq_folder`:The absolute path to the folder in which all raw fastq files are stored.
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
- `dedup_umi`: Whether or not the dedupliucation of UMIs step should takingTrue or False
- `bc_pattern`: If umi deduplication is taking place, which barcode pattern should be given to UMI tools

# Commands

 Pipeline commands (must be run alongside samples command)

- `all`: runs all of the pipelines
- `align`: merges (if indicated in parameters.csv) and aligns reads to ref genome
- `identify`: identifies nominated off-target sites
- `visualize`:annotates identified sites and creates alignment .svg image
- `coverage`: create .bam file, alignment histogram and stats file for the identified sites
- `variants`: *currently* updating

Replicates can be combined seperate calling `replicate_combiner.py`, see below

### Samples

- `all`: runs all samples in manifest
- `sample_name`: runs the specific sample indicated

## Stand alone scripts

### Replicate Combiner

Normalizes sample replicates and combines counts into one .csv file. Also creates venn diagram, scatter plot and swarm plot

```
python path-to/IGU_CHANGEseq/changeseq/replicate_combiner.py --sample1 test_replicate1 --sample2 test_replicate2 --name test --output path-to/TEST_DATA/  --read_threshold 6
```

### Multisample combined

in progess

# Pipeline Output
When running the full pipeline, the results of each step are outputted to the `output_folder` in a separate folder for each step. The output folders and their respective contents are as follows:

.....
