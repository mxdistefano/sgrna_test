import pandas as pd
import sys

long_data=sys.argv[1]

colnames = ['gene_name','hgnc_id','tpm_unstranded','sample_name']
out_matrix=pd.read_csv(long_data,sep="\t",names=colnames)
out_matrix.pivot_table(index=["gene_name", "hgnc_id"], 
                columns='sample_name', 
                values='tpm_unstranded').to_csv("TPM_expression_matrix_TCGA-BRCA_selected_samples.tsv",sep="\t")