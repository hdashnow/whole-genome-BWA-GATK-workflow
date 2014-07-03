load 'pipeline_stages_config.groovy'

run {
    "%.fastq.gz" * [ fastqc ] + 
    "%_*.fastq.gz" * [ 
        align_bwa + index_bam +
        dedup + index_bam + 
        realignIntervals + realign + index_bam +
        recal_count + recal + index_bam +
        call_variants
        ]
}
