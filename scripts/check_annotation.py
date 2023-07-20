"""
this script takes as argument the mapping coordinates annotated
and check in what sequences the gene name predicted matches the
truth gene encode in the sequence id  
input : arg1: output of bedtools intersect annotation (tsv file)  
		arg2: hgnc_id table from HUGO

output : -matching_entries.csv tsv file containing well annotated sequences
		 -discordant_entries.csv tsv file containing sequences whose gene 
		 annotation did not match its true gene

example:
python check_annotation.py mapping_annotated.bed
"""

import sys
import pandas as pd

annotated_file=sys.argv[1]
hgnc_table_file=sys.argv[2]

"""
parse gene name from gencode annotation.
Gene name position is different for each category
"""


def parse_gene_name(gencode_ann):
    gencode_ann_splitted=gencode_ann.split(";")
    for block in gencode_ann_splitted:
        if "gene_name" in block:
            return block.split(" ")[2].replace('"','')
    return None

def parse_hgnc_id(gencode_ann):
    gencode_ann_splitted=gencode_ann.split(";")
    for block in gencode_ann_splitted:
        if "hgnc_id" in block:
            return block.split(" ")[2].replace('"','')
    #if not in annotation description search in dict
    return rev_hgnc_id_dict.get(parse_gene_name(gencode_ann),None)

#build dict to translate Gene symbol to hgnc id
def build_hgnc_id_dict(hgnc_table):
	##join all gene name aliasses 
	hgnc_values=hgnc_table[['symbol','alias_symbol','alias_name','prev_symbol']].apply(
    lambda x: '|'.join(x.dropna().astype(str)),
    axis=1)
	
	hgnc_id_dict = dict()
	hgnc_ids=hgnc_table["hgnc_id"]

	for i in range(len(hgnc_values)):
	    for v in hgnc_values[i].replace(" ", "").split("|"):
	        if hgnc_ids[i] in hgnc_id_dict:
	            hgnc_id_dict[hgnc_ids[i]].append(v)
	        else:
	            hgnc_id_dict[hgnc_ids[i]] = [v]

	#reverse dict for mapping           
	rev_hgnc_id_dict = {}
	for k,v in hgnc_id_dict.items():
	    for i in v:
	        rev_hgnc_id_dict[i] = k
	return rev_hgnc_id_dict

if __name__ == "__main__":
	hgnc_table=pd.read_csv(hgnc_table_file,sep="\t", usecols =['symbol','hgnc_id','alias_symbol','alias_name','prev_symbol','ensembl_gene_id'])

	rev_hgnc_id_dict=build_hgnc_id_dict(hgnc_table)

	#annotated_file='full_library_mapped_bestmatch_formated_ann.bed'
	mapping=pd.read_csv(annotated_file,sep="\t",header=None)
	#rename some columns and drop unwanted
	mapping=mapping.rename(columns = {0:'chr',1:'start',2:'end',3:'id',6:'region',10:'strand',12:'annotation'})
	col2drop=[x for x in mapping.columns if str(x).isdigit()]
	mapping.drop(mapping.columns[col2drop],axis=1,inplace=True)

	#parse gene name from fasta id/description
	mapping["true_gene"]=mapping['id'].str.split("|",expand=True)[2].str.split("_",expand=True)[0]
	#parse gene name from gff annotation
	mapping["predicted_gene"]=mapping['annotation'].apply(lambda x : parse_gene_name(x))
	#map gene name to get hgnc_id
	mapping['hgnc_id']=mapping['true_gene'].map(rev_hgnc_id_dict).fillna("NA")
	mapping["hgnc_id_predicted"]=mapping['annotation'].apply(lambda x : parse_hgnc_id(x))

	#remove duplicate hits
	mapping=mapping.groupby(['id','predicted_gene']).first().reset_index()
	#get sequences with matching genes from fasta description and mapping
	matching=mapping[mapping.apply(lambda x: x["hgnc_id_predicted"] == x["hgnc_id"], axis=1)].reset_index()
	#select discordant sequences
	no_matching=mapping[mapping.apply(lambda x: x["hgnc_id_predicted"] != x["hgnc_id"], axis=1)].reset_index()
	#filter secondary hits of matching sequences that match elsewhere
	unique_no_matching=no_matching.loc[~no_matching['id'].isin(matching['id']), ]

	matching.to_csv("matching_sequences.csv",sep="\t")
	unique_no_matching.to_csv("no_matching_sequences.csv",sep="\t")
