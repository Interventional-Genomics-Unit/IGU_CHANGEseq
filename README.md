# IGI CHANGE-seq & CHANGE-seq BE
 An extension of CHANGE-seq for Off-Target Site Nomination
 https://github.com/tsailabSJ/changeseq
 
Added Features:

* updated to Python3
* ability to limit search with bulges
* ability to limit search with edit distance
* simple .csv sample manifest
* removal of TN5 ME adapters and illumina adapter (unmerged version only)
* replicate combination
* normalization of replicates - (median-of-ratios or rpm)
* QC metrics of alignments and fastq circularization read composition
* identified site bam file for use with IGV
* genomic annotation to sites
* more plotting options - venn diagram, swarmplot and more
* Fixed variant analysis bug
* coming soon: dual BE and Cas pipline 


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
cd /your-path/IGU_CHANGEseq/
python changeseq/changeseq.py makefiles
```
 
## Usage

The Tsai lab change-seq pipeline requires a manifest yaml file specifying input files, output directory, and pipeline parameters. 
In the IGU version this file is split into t two csv files 1) a sample manifest file <code>--manifest</code> 2) an optional settings 
file <code>--settings</code>
Once the csv files are created, users can simply run the changeseq.py with the parameters indicating the 
raw fastq directory <code>--raw_fastq_folder</code> and main analysis directory <code>--analysis_folder</code>


The bash wrapper is initiated with the parameters.csv, manifest.csv, command, and sample_name

```
cd /your-path/IGU_CHANGEseq/
python changeseq/changeseq.py all --analysis_folder /your-path/IGU_CHANGEseq/changeseq/test/Standard_Output --raw_fastq_folder /your-path/IGU_CHANGEseq/changeseq/test/input --settings unmerged_parameters.csv --manifest manifest.csv --sample all
```

### Writing the Sample Manifest
The original manifest.yaml is no longer used in this version. The prior manifest.yaml is now split into two .csv input files. 
The manifest file has the following .csv file 


- `sequencing_sample_name`: sample name as indicated in the raw fastq (Leave out the ("_R1_001")
- `sample_name`: The sample name as you wish to see on the output file names  
- `control_sequencing_sample_name`: No RNP control sample name as indicated in the raw fastq (Leave out the ("_R1_001")
- `target`: Target sequence for that sample. Accepts degenerate bases.
- `description`: A brief description of the sample. Usually the cell or gDNA Source name
- `replicate_group_name`: (optional) an identifier that designates replicate samples. This must be present to run normalization  
     

 ### Parameters and Setting
 This file is a .csv with the first column the parameter name and the second column the parameter name. 
 The parameters are now stored internally and only need to be given if different than default
 See IGU_CHANGEseq/changeseq/test/parameters.csv


- `bwa`: The absolute path to the `bwa` executable
- `samtools`: The absolute path to the `samtools` executable
- `cutadapt`: The absolute path to the `cutadapt` executable
- `merged_analysis`: Whether or not the paired read merging step should taking (highly recommend not using) default = False
- `read_threshold`: The minimum number of reads at a location for that location to be called as a site. We recommend leaving it to the default value of 4.
- `window_size`: Size of the sliding window, we recommend leaving it to the default value of 30.
- `mapq_threshold`: Minimum read mapping quality score. We recommend leaving it to the default value of 40.
- `start_threshold`: Tolerance for breakpoint location. We recommend leaving it to the default value of 1.
- `gap_threshold`: Distance between breakpoints. We recommend leaving it to the default value of 3 for Cas9.
- `mismatch_threshold`: Number of tolerated mismatches in the fuzzy target search step. We recommend leaving it to the default value of 6.
- `edist_threshold`: Maximum edit (levenshtein) distance in the fuzzy target search step. default=7
- `bulge_threshold`: Number of tolerated indels (bulges) in the fuzzy target search step. default=1
- `read_length`: Fastq file read length, default is 151.
- `PAM`: PAM sequence, default is NGG.
- `search_radius`: in addition to the window size, how far away from the cut site to search for search
- `normalize`: normalization method to use for replicates "none", "median" or "rpm", default=median

For BE analysis
- `BEsearch_radius`: In addition to the window size, how far away from the overlap to search for search 
- `BEmodel_min_overlap`: Min. overlap between reads required. Ex. a value of 2 means the nCas9 would nick 2bp from the deaminase. Highly recommended keeping this at 2.
- `BEmodel_max_overlap`: Max overlap between reads required. Ex. a value of 15 means the nCas9 would nick 15bp from the deaminase. 


### Pipeline commands (must be run alongside samples command)

- `all`: runs all of the pipelines
- `align`: merges (if indicated in parameters.csv), trims Tn5 adapters and aligns reads to ref genome
- `identify`: identifies nominated off-target sites and annotates sites
- `visualize`:annotates identified sites and creates alignment .svg image
- `qc`: create .bam file for identified sites for IGV, alignment histogram and alignment stats file and additional QC about fastq 
- `analyze`: joines replicates if given, normalizes reads and creates additional plots
- `variants`: *currently* updating

### Parameters

- `analysis_folder`: The absolute path to the folder in which all pipeline outputs will be saved.
- `raw_fastq_folder`:The absolute path to the folder in which all raw fastq files are stored.
- `all`: runs all samples in manifest
- `sample_name`: runs the specific sample indicated
- `base-editing`: run CHANGE-seq BE (T or F). default=F


## Pipeline Results
When running the full pipeline, the results of each step are outputted to the `output_folder` in a separate folder for each step. The output folders and their respective contents are as follows:

```
project_directory/
│
├── raw_fqs/                        # User created input folder of FASTQs
│
├──aligned/                         # aligned reads
│
├──preprocessed/                    # Cutadapt & fastp output for FASTQ trimming 
│
├── raw_output/                     # Pre-normalized raw read counts
│   └── visualization/              # Alignment plots of raw read counts
│
├── qc/
│   ├── alignment_files/            # BAM stats and a subset of identified sites (great for IGV!)
│   └── QC_reports/                 # Reports: % trimmed, % aligned, read counts, etc.
│
└── processed_output/               # Processed + normalized reads
    ├── joined_tables/              # Detailed replicate-level summary tables
    └── visualizations/             # Swarm plots, scatter plots, alignment plots, etc.

```