nextflow.enable.dsl=2
hgnc_id_file=file(params.hgnc_table)
ref_file=file(params.reference)

process mapping {
	publishDir "$params.publish_dir", mode: 'copy'
	
	input:
	path input_fasta
	//path ref

	output:
	path 'fasta_mapped.psl' , emit: mapped_fasta

	script:
	"""	
	pblat -threads=$params.mapping_threads stepSize=5 tileSize=11 -fine -noHead -minScore=19 -minIdentity=95 ${ref_file} ${input_fasta} fasta_mapped.psl
	"""
}

process format_mapping_output {
	input:
	path mapping_out

	output:
	path 'fasta_mapped_formatted.psl', emit: mapped_fasta_formatted

	script:
	"""
	cut -d'	' -f10,14,16,17  ${mapping_out} | awk '{ print \$2 "\t" \$3 "\t" \$4 "\t" \$1}' > fasta_mapped_formatted.psl
	"""
}

process annotate {
	publishDir "$params.publish_dir", mode: 'copy'

	input:
	path fasta_mapped
	path gencode_gtf

	output:
	path 'fasta_mapped_annotated.bed', emit: annotated_file
	
	script:
	"""
	bedtools.static intersect -a ${fasta_mapped} -b ${gencode_gtf} -wao > fasta_mapped_annotated.bed
	"""
}

process check_annotation{
	publishDir "$params.publish_dir", mode: 'copy'

	input:
	path check_ann_script
	path annotated_file

	output:
	path 'matching_sequences.csv', emit: matching_seqs
	path 'no_matching_sequences.csv', emit: no_matching_seqs

	script:
	"""
	python3 ${check_ann_script} ${annotated_file} ${hgnc_id_file}
	"""
}


process download_tcga_data{
	publishDir "$params.publish_dir", mode: 'copy'

	input: 
	tuple val(sample_id), val(uuid)

	output:
    path '*.tsv'

	script:
	"""
	gdc-client download --no-related-files $uuid
	mv ${uuid}/*.tsv ${sample_id}_${uuid}.tsv
	"""
}

process filter_tcga_data{
	publishDir "${params.publish_dir}", mode: 'copy'
	input: 
	path filter_tcga_script
	path tcga_file
	path matching_seqs

	output:
	path '*.tsv'

	script:
	"""
	python3 ${filter_tcga_script} ${matching_seqs} ${tcga_file} ${hgnc_id_file} 
	"""
}

process format_final_data{
	publishDir "${params.publish_dir}", mode: 'copy'
	input: 
	path long_to_wide_matrix_script
	path long_data

	output:
	path '*.tsv'

	script:
	"""
	python3 ${long_to_wide_matrix_script} ${long_data}
	"""
}


workflow {
	input_fasta_file=file(params.input_fasta)
	gencode_gtf_file=file(params.gtf_file)
	check_ann_script=file(params.check_ann_script)
	filter_tcga_script=file(params.filter_tcga_script)
	long_to_wide_matrix_script=file(params.long_to_wide_matrix)

	Channel.from(params.tcga_samples).set{samples_to_download}
	//sample_uuids=samples_to_download.map { it[1] }

    mapping(input_fasta_file)
    format_mapping_output(mapping.out.mapped_fasta)
    annotate(format_mapping_output.out.mapped_fasta_formatted,gencode_gtf_file)
    check_annotation(check_ann_script,annotate.out.annotated_file)
	tcga_files=download_tcga_data(samples_to_download)
    filtered_data=filter_tcga_data(filter_tcga_script,tcga_files,check_annotation.out.matching_seqs)
    //collect output files to make one expression matrix with both samples as columns
    filtered_data.collectFile(name: 'all.tsv',storeDir: params.publish_dir).set{concat_filtered_data}
    //convert to wide format
    format_final_data(long_to_wide_matrix_script,concat_filtered_data)
}
