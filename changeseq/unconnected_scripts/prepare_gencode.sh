#bash

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz -O genecode_annotation.gtf.gz
gunzip genecode_annotation.gtf.gz
awk '$2 == "HAVANA" { print }' genecode_annotation.gtf > genecode_annotation_processed.gtf
sed -i 's/;//g' genecode_annotation_processed.gtf
sed -i 's/"//g' genecode_annotation_processed.gtf

#cat genecode_annotation_hg38.gtf | awk '{OFS="\t"} {print $1,$4-1,$5,$3,$10,$16,$7}' | tr -d '";' > genecode_hg38.bed
awk '$3 == "gene"&& $12 == "protein_coding" { print $1 "\t" $4 "\t" $5 "\t" $3 "\t" $10 "\t" $14}' genecode_annotation_processed.gtf > GENE_genecode_annotation.bed
awk '$3 == "exon"&& $14 == "protein_coding" { print $1 "\t" $4 "\t" $5 "\t" $3}' genecode_annotation_processed.gtf > EXON_genecode_annotation.bed
awk '$3 == "UTR"&& $14 == "protein_coding" { print $1 "\t" $4 "\t" $5 "\t" $3}' genecode_annotation_processed.gtf > UTR_genecode_annotation.bed
awk '$3 == "start_codon"&& $14 == "protein_coding" { print $1 "\t" $4 "\t" $5 "\t" $3}' genecode_annotation_processed.gtf > START_CODON_genecode_annotation.bed
awk '$3 == "stop_codon"&& $14 == "protein_coding" { print $1 "\t" $4 "\t" $5 "\t" $3}' genecode_annotation_processed.gtf > STOP_CODON_genecode_annotation.bed

awk '$3 == "transcript" && $2 == "HAVANA" { print }' genecode_annotation.gtf > TRANSCRIPT_genecode_annotation.gtf
awk '$3 == "CDS" && $2 == "HAVANA" { print }' genecode_annotation.gtf > CDS_genecode_annotation.gtf


#less genecode_annotation_processed.gtf