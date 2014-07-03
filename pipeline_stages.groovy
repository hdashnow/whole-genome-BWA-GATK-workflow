////////////////////////////////////////
BASE="/vlsci/VR0320/shared/production"
// Sequencing platform
PLATFORM="illumina"
// Where you installed GATK
GATK="gatk"
// Set location of you reference files 
REFBASE="hg19"
REF="$REFBASE/gatk.ucsc.hg19.fasta"
DBSNP="$REFBASE/dbsnp_132.hg19.vcf"
GOLD_STANDARD_INDELS="$REFBASE/Mills_and_1000G_gold_standard.indels.b37.chr.vcf"
////////////////////////////////////////


fastqc = {
    doc "Run FASTQC on raw reads"
    output.dir = "fastqc"
    transform('.fastq.gz')  to('_fastqc.zip') {
        exec "fastqc -o ${output.dir} $inputs.gz"
    }
}

align_bwa = { 
    transform("bam","sai","sai") {
        output.dir = "align"        
        
        def SAMPLE = "NA10830"
        def ID = "1"

        msg "Aligning $input ..."

        multi "gunzip -c $input1 | bwa aln $HGFA - > $output2",
              "gunzip -c $input2 | bwa aln $HGFA - > $output3"
        
        msg "Running bwa sampe"
        exec """
            bwa sampe -r "@RG\\tID:${ID}\\tPL:$PLATFORM\\tPU:None\\tLB:None\\tSM:${SAMPLE}" -n 5 -N 20 $HGFA $output2  $output3 $input1 $input2 | samtools view -bSu - | samtools sort - $output1.prefix
        """
   }
}

index_bam = {

    doc "Create BAM file index"
    
    //make sure index is in same folder as bam
    output.dir=file(input.bam).absoluteFile.parentFile.absolutePath 
    transform("bam") to ("bam.bai") {
        exec "samtools index $input.bam","index_bam"
    }
    forward input
}

realignIntervals = {
    doc "Select regions for realignment with GATK in 'realign' stage"
    output.dir="align"
    exec """
        java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
            -T RealignerTargetCreator 
            -R $REF 
            -I $input.bam 
            --known $GOLD_STANDARD_INDELS 
            -o $output.intervals
    """, "realign_target_creator"
}

realign = {
    doc "Apply GATK local realignment to indevals selected in stage 'realignIntervals' "
    output.dir="align"
    exec """
        java -Xmx5g -jar $GATK/GenomeAnalysisTK.jar 
             -T IndelRealigner 
             -R $REF 
             -I $input.bam 
             -targetIntervals $input.intervals 
             -o $output.bam
    ""","local_realign"
}

dedup = {
    doc "Remove PCR duplicate reads from BAM"
    output.dir="align"
    exec """
        MarkDuplicates
             INPUT=$input.bam 
             REMOVE_DUPLICATES=true 
             VALIDATION_STRINGENCY=LENIENT 
             AS=true 
             METRICS_FILE=$output.metrics
             OUTPUT=$output.bam
    ""","MarkDuplicates"
}

recal_count = {
    doc "Recalibrate base qualities in a BAM file using observed error rates"
    output.dir="align"
    INDEL_QUALS=""
    // To use lite version of GATK uncomment below
    // INDEL_QUALS="--disable_indel_quals"

    exec """
        java -Xmx5g -jar $GATK/GenomeAnalysisTK.jar 
             -T BaseRecalibrator 
             -I $input.bam 
             -R $REF 
             --knownSites $DBSNP $INDEL_QUALS
             -l INFO 
             -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate 
             -o $output.counts
    """, "recalibrate_bam"
}

recal = {
    doc "Apply recalibration quality adjustments"
    output.dir="align"
    exec """
        java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
           -T PrintReads 
           -I $input.bam 
           -BQSR $input.counts 
           -R $REF 
           -l INFO 
           -o $output.bam
        """, "recalibrate_bam"
}

call_variants = {
    doc "Call SNPs/SNVs using GATK Unified Genotyper"
    output.dir="variants"

    // Default values from Broad
    var call_conf:5.0, 
        emit_conf:5.0

    transform("bam","bam") to("metrics","vcf") {
        exec """
            java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar 
                   -T UnifiedGenotyper 
                   -R $REF 
                   -I ${inputs.bam.withFlag("-I")} 
                   -nt 4
                   --dbsnp $DBSNP 
                   -stand_call_conf $call_conf -stand_emit_conf $emit_conf
                   -dcov 1600 
                   -l INFO 
                   -A AlleleBalance -A Coverage -A FisherStrand 
                   -glm BOTH
                   -metrics $output.metrics
                   -o $output.vcf
            ""","gatk_call_variants"
    }
}
