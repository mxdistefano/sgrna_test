# sgrna_test

Pipeline for processing sgRNA sequences as indicated. Consist in the following steps:
1-Mapping: input fasta file to GRCh38 genome
	-software used: pblat . Only allowed one mismatch. 
	-Download GRCh38 fasta reference from https://www.gencodegenes.org/human/ and reference path in configuration file 'nextflow.config'

2-Annotation: using a gtf associated to the chosen reference genome:
	-software used: bedtools
	-Download gtf file from https://www.gencodegenes.org/human/ and reference path in configuration file 'nextflow.config'

3-check_annotation: for each sequence compare the gene contained in fasta description and the gene obtained in the annotation. 
	-software used: a python script was developed to  translate each gene name to its corresponding hgnc_id because some sgRNA sequences had non canonical names and make the corresponding data manipulation and comparisons.

This process generates two outputs:
	-matching_sequences.csv : sequences correctly annotated. This file will be passed to the next process
	-no_matching_sequences.csv : sequences with discordant gene annotation. Should be manually inspected.
Both files contain the genomic coordinates of the mapping and the annotated gene.

4-download_tcga_data: process to download the requested TCGA-BRCA data. The sample ID and uuid of the experiment mus be added in 'nextflow.config' to make the download possible.
	software used: gdc-client

5-filter_tcga_data: Filter the previously downloaded files to obtain only the genes from the “matching sequences”
	-software used: a python script was developed to filter the files and format an expression matrix in tpm as output. 
Final output: TPM_expression_matrix_TCGA.tsv

Requirement:
All the required software has been packed in the following docker image: 
		https://hub.docker.com/r/mxdistefano/sgrna_tools/tags  (docker pull mxdistefano/sgrna_tools:latest)
		or it can be build from the dockerfile in the folder 'docker_image'

Reference sequence and annotation can be obtained from:	https://www.gencodegenes.org/human/
