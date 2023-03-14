
process combine_fastqs {
	input:
	val sample_id
	path read_dir

	output:
	tuple val(sample_id), path("${sample_id}_1.fastq.gz"), path("${sample_id}_2.fastq.gz")
	
	shell:
	'''
	echo !{sample_id}
	forward_files="\$(find !{read_dir}/ | grep -E "!{sample_id}_.*(_1)(\\.fastq\\.gz|\\.fq\\.gz)" | sort)"
	reverse_files="\$(find !{read_dir}/ | grep -E "!{sample_id}_.*(_2)(\\.fastq\\.gz|\\.fq\\.gz)" | sort)"
	if [[ "\$(printf '%s\\n' "${forward_files[@]}" | wc -l)" -gt 1 ]]
	then
		unpigz -p 8 -c \${forward_files[@]} > "!{sample_id}_1.fastq"
		unpigz -p 8 -c \${reverse_files[@]} > "!{sample_id}_2.fastq"
		pigz -p 8 !{sample_id}_1.fastq
		pigz -p 8 !{sample_id}_2.fastq
	else
		ln -s ${forward_files} "!{sample_id}_1.fastq.gz"
		ln -s ${reverse_files} "!{sample_id}_2.fastq.gz"
	fi
	'''
}


