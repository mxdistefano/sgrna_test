import os, sys, glob
import pandas as pd	
from check_annotation import build_hgnc_id_dict
if __name__ == "__main__":
    matching_reads_file=sys.argv[1]
    uuid_tcga_file=sys.argv[2]
    hgnc_table_file=sys.argv[3]
    #out_file_name=sys.argv[4]

    counts_metric='tpm_unstranded'
    
    hgnc_table=pd.read_csv(hgnc_table_file,sep="\t", usecols =['symbol','hgnc_id','alias_symbol','alias_name','prev_symbol','ensembl_gene_id'])

    hgnc_id_dict= build_hgnc_id_dict(hgnc_table)

    #select unique list of genes from sequences that matched the annotation
    matching_reads=pd.read_csv(matching_reads_file,sep="\t")
    unique_gene_set=matching_reads['hgnc_id'].unique()
    unique_gene_set_name=matching_reads['true_gene'].unique()

    #Selected tpm_unstranded column, buit it could be changed according to interest
    cols_of_interest=['hgnc_id','gene_name',counts_metric]

    uuid=uuid_tcga_file.split("_")[0]
    expr_matrix=pd.read_csv(uuid_tcga_file,sep="\t",comment="#")
    expr_matrix["hgnc_id"]=expr_matrix['gene_name'].map(hgnc_id_dict).fillna("NA")
    #search for genes by hgnc_id. If not found search by gene name
    expr_matrix_matching=expr_matrix[expr_matrix.apply(lambda x: True if x['hgnc_id'] in unique_gene_set else( True if x['gene_name'] in unique_gene_set_name else False)
                                                                    ,axis=1)]

    #if outfile exist, other uuid was processed before. Then open file to append column for this one.
    #else create an empty one 
    index_cols=["gene_name", "hgnc_id"]
    expr_matrix_matching=expr_matrix_matching[cols_of_interest].set_index(index_cols)
    expr_matrix_matching=expr_matrix_matching.rename(columns = {counts_metric:uuid+'_'+counts_metric})
    expr_matrix_matching["sample"]=uuid
    #if os.path.isfile(out_file_name): 
    #    out_matrix=pd.read_csv(out_file_name,usecols=cols_of_interest,comment="#",index_col=index_cols)
    #    df2 = pd.merge(out_matrix,expr_matrix_matching , left_index=True, right_index=True)
    #else:
    #    out_matrix=expr_matrix_matching
    #out_matrix.to_csv(uuid+_+"filtered_genes.tsv",sep="\t")
    expr_matrix_matching.to_csv(uuid+'_filtered_genes.tsv',sep="\t",header=False)