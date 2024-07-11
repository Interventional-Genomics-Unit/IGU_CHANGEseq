# IGU CHANGE-seq
 An extension of CHANGE-seq for Off-Target Site Nomination 
 
 

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
 
