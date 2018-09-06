#!/usr/bin/env nextflow

/*------------------------------------------------------------------
stp-bioinformatics@crick.ac.uk Nextflow script for acient DNASeq analysis
------------------------------------------------------------------*/

// print command-line help
def print_help( ) {
  log.info 'aDNASeq-byBABS ~ ancient DNA-Seq Analysis Pipeline'
  log.info '--------------------------------------------------'
  log.info 'Run an aDNASeq analyse on a set of fastq files.'
  log.info ''
  log.info 'Usage using wrapper script: '
  log.info '    aDNASeq-byBABS --outdir [ output directory path ] --genome [ hg38 ] --design [ path to design file ] -profile [ profile name ]'
  log.info ''
  log.info 'Options:'
  log.info '    --help                      Show this message and exit.'
  log.info ''
  log.info 'Parameters:'
  log.info '    --genome                    Name of genome bundle defined in config.'
  log.info ''
  log.info '    --design CSV_FILE           Comma separted file containing information'
  log.info '                                about the samples (see README).'
  log.info ''
  log.info '    -profile                    Compute environment profile to use for processing'
  log.info '                                [babs|conda|none] (see nextflow.config)'
  log.info ''
}

//print help
if (params.help) {
  print_help( )
  exit( 0 )
}

// CHECK THAT NEXTFLOW VERSION IS UP TO DATE ENOUGH
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "============================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

// the java modules
import java.nio.file.Paths
import java.nio.file.Files

// exotic source modules
import org.yaml.snakeyaml.Yaml
import org.codehaus.groovy.runtime.StringGroovyMethods

version = 1.0

/* verify parameters and paths
 */

// check genome name and file paths exist and set in params
if( params.genome ) {
  if( ! (genome_paths = params.genomes[ params.genome ]) ) {
    throw new Exception( params.genome + " is an invalid genome name, must be one of " + params.genomes )
  }else{
    for( String key : genome_paths.keySet() ){
      if( !new File( genome_paths.get(key) ).exists() ) {
	throw new Exception( key + " reference genome file does NOT exist.\n" + genome_paths.get( key ) )
      }else{
	params.put( key, genome_paths.get(key) )
      }
    }
  }
}else{
  println( "Error: --genome not set\n" )
  print_help()
  exit( 1 )
}

// check design file exists
if( params.design ) {
  if( !(f = new File( params.design ) ).exists() ) {
    throw new Exception( "design file " + f.getPath() + " NOT found" )
  }
} else {
  println( "Error: --design not set\n" )
  print_help()
  exit( 1 )
}

/* Set output directory to current working directory.
   It must be like this so all nextflow reports are written to 
   the same place as the analysis results.
*/
params.outdir = "./"

// set multiqc config path
params.multiqc_config = new File( "${workflow.projectDir}/conf/multiqc_config.yml" )
if( !params.multiqc_config.exists( ) ) {
  throw new Exception( "multiqc config file " + params.multiqc_config + " NOT found" )
}

/*
initialise primary channel from design file
*/
Channel
  .fromPath( params.design )
  .splitCsv( header:true )
  .map { row ->
    [ row.library,
      row.individual,
      file( row.R1 ),
      file( row.R2 ) ] }
  .set{ fastq_ch } 

/* mergeFastq.
   Merge fastqs from the same library
*/
fastq_ch
    .groupTuple()
    .set{ fastq_mergeFastq_ch }

process merge_fastq {
    
  tag{ library }
    
  input: set library, indiv, fastqR1, fastqR2 from fastq_mergeFastq_ch

  output: set library, indiv, "*fastq.gz" into mergedFastq_ch

  script:
  def fastqR1String = fastqR1.join(" ")
  def fastqR2String = fastqR2.join(" ")
  def fastqcount = fastqR1.size()
  """
  if [ ${fastqcount} == 1 ]
  then
    ln -s ${fastqR1String} ${library}_R1.fastq.gz
    ln -s ${fastqR2String} ${library}_R2.fastq.gz
  else
    zcat ${fastqR1String} | gzip > ${library}_R1.fastq.gz
    zcat ${fastqR2String} | gzip > ${library}_R2.fastq.gz
  fi
  """
}

/* SeqPrep: 
    Remove adapter and poor quality sequence from the 3' & 5' end of reads.
    The maximum read error rate is set to 10% (-e 0.1), 
    the quality cutoff is set at 10 (-q 10),
    the minimum acceptable read length after trimming is 20 (-m 20)
    and the minimum adapter overlap is 1 (-O 1).
*/

process seqprep {
  
  tag { library }

  cache 'deep'
  
  publishDir params.outdir, mode: 'copy', overwrite: true,
  saveAs: {filename ->
    if( filename.endsWith(".log")) "log/$filename"
    else if( filename.endsWith(".stats")) "qc/$filename" }
  
  input: set library, indiv, fastqfiles from mergedFastq_ch
  
  output: set library, indiv, "*SeqPrep.fastq.gz" into seqprep_ch
  output: file "*.log" into seqprep_log_ch
  output: file "*.stats" into seqprep_stats_ch
  
  script:
  """
  SeqPrep \
    -f ${fastqfiles[0]} \
    -r ${fastqfiles[1]} \
    -1 ${library}.SeqPrep_R1.fastq.gz \
    -2 ${library}_SeqPrep_R2.fastq.gz \
    -A CTGTCTCTTATA -B CTGTCTCTTATA \
    -s ${library}.SeqPrep.fastq.gz > ${library}.SeqPrep.log 2>&1
  cp ${library}.SeqPrep.log ${library}.SeqPrep.stats
  """
}

/* bwa aln
 */
process bwa_aln {
  
  tag{ library }
  
  publishDir params.outdir, mode: 'copy', overwrite: true,
  saveAs: {filename ->
    if( filename.endsWith(".log")) "log/$filename" }
	
  input:
    set library, indiv, fastq from seqprep_ch
    //file index from bwaIndexCh.collect( )
  
  output:
	set library, indiv, fastq, "*.sai" into bwa_aln_ch
	file "*.log" into bwa_aln_log_ch
    
  script:
  """
  bwa \
    aln -t ${task.cpus} \
    -l 16500 \
    -n 0.01 \
    -o 2 \
    ${params.bwa_idx} \
    ${fastq} > ${library}.sai 2> ${library}_bwa_aln_.log
  """
}

process bwa_bam {

    tag{ library }

    input:
	set library, indiv, fastq, sais from bwa_aln_ch
    //file index from bwaIndexCh.collect( )
  
    output:
    set library, indiv, "*.bam" into bwa_bam_markdups_ch, bwa_bam_qc_ch
    file "*.log" into bwa_bam_log_ch
    
    script:

    def rg="\'@RG\\tID:${library}\\tSM:${library}\\tPL:illumina\\tLB:1\\tPU:1\'"
    """
    bwa samse \
    -r $rg ${params.bwa_idx} \
    ${sais} \
    ${fastq} > ${library}.sam 2> ${library}_bwa_sam.log
    samtools view -b -o ${library}_unsorted.bam ${library}.sam
    samtools sort -o ${library}.bam ${library}_unsorted.bam
    samtools index ${library}.bam
    rm ${library}.sam
    rm ${library}_unsorted.bam
    """
}

/* Picard Mark Duplicates
 */
process remove_duplicates {
  
  tag{ library }

  label 'picard'
  
  publishDir params.outdir, mode: 'copy', overwrite: true,
  saveAs: {filename ->
	if( filename.endsWith(".bam")) "results/alignments/$filename"
	else if (filename.endsWith(".bai")) "results/alignments/$filename"
	else if (filename.endsWith("log")) "pipeline/logs/$filename"
	else if (filename.endsWith("MarkDuplicates.metrics.txt")) "qc/data/$filename"
	else if (filename.endsWith("flagstat" )) "qc/data/$filename"
	else if (filename.endsWith("idxstats" )) "qc/data/$filename" }

  input:
  set library, indiv, bamfile from bwa_bam_markdups_ch

  output:
  set library, indiv, "*.bam" into(
    remove_duplicates_ch,
    remove_duplicates_pmdtools_metrics_lib_ch,
    remove_duplicates_oxoG_metrics_lib_ch,
    remove_duplicates_wgs_metrics_lib_ch,
    remove_duplicates_insert_size_lib_ch,
    remove_duplicates_collect_alignment_metrics_library_ch,
    remove_duplicates_bam_to_vcf_library_ch )
  file "*.log" into remove_duplicates_log_ch
  file "*.MarkDuplicates.metrics.txt" into remove_duplicates_metrics_ch
  file "*.flagstat" into remove_duplicates_flagstat_ch
  file "*.idxstats" into remove_duplicates_idxstats_ch

  script:
  def java_mem = new String( task.memory.toString() ).replace( " GB", "" )
  """
  java -Xmx${java_mem}g -Djava.io.tmpdir=./ \
    -jar \${EBROOTPICARD}/picard.jar MarkDuplicates \
    VALIDATION_STRINGENCY=LENIENT \
    REMOVE_DUPLICATES=true \
    ASSUME_SORTED=true \
    TMP_DIR=./ \
    INPUT=${bamfile} \
    OUTPUT=${library}.dedup.bam \
    METRICS_FILE=${library}.MarkDuplicates.metrics.txt \
    >> ${library}.MarkDuplicates.log 2>&1
  samtools index ${library}.dedup.bam
  samtools flagstat ${library}.dedup.bam > ${library}.bam.flagstat
  samtools idxstats ${library}.dedup.bam > ${library}.bam.idxstats
  """
}

process pmdtools_lib {

  tag{ library }

  errorStrategy 'ignore'

  publishDir params.outdir, mode: 'copy', overwrite: true,
  saveAs: {filename ->
    if (filename.endsWith("log")) "pipeline/logs/$filename"
    else if (filename.endsWith("stats")) "qc/data/$filename" }
	
  input:
  set library, indiv, bamfile from remove_duplicates_pmdtools_metrics_lib_ch

  output:
  file "*.log" into pmd_metrics_lib_log_ch
  file "*.stats" into pmd_metrics_lib_stats_ch
  
  script:
  """
  samtools view -q30 ${bamfile} | pmdtools --first \
  2> ${library}.pmdtools_first.log
  cp ${library}.pmdtools_first.log ${library}.pmdtools_first.stats
  ## To compute deamination-derived damage patterns separating CpG and non-CpG sites
  samtools view -q30 ${bamfile} | pmdtools --platypus > ${library}.pmdtools_platypus.stats \
  2>  ${library}.pmdtools_platypus.log
  samtools view -q30 ${bamfile} | pmdtools --deamination --range 30 --CpG \
  > ${library}.pmdtools_deamination.stats 2> ${library}.pmdtools_deamination.log
  """
}

process oxoG_metrics_lib {

  tag{ library }

  label 'picard'

  publishDir params.outdir, mode: 'copy', overwrite: true,
  saveAs: {filename ->
    if (filename.endsWith("log")) "pipeline/logs/$filename"
    else if (filename.endsWith("stats")) "qc/data/$filename" }
	
  input:
  set library, indiv, bamfile from remove_duplicates_oxoG_metrics_lib_ch

  output:
  file "*.log" into oxoG_metrics_lib_log_ch
  file "*.stats" into oxoG_metrics_lib_stats_ch
  
  script:
  def java_mem = new String( task.memory.toString() ).replace( " GB", "" )
  """
  java -Xmx${java_mem}g -Djava.io.tmpdir=./ \
    -jar \${EBROOTPICARD}/picard.jar CollectOxoGMetrics \
    R=${params.fasta} \
    I=${bamfile} \
    O=${library}.oxoG_metrics_lib.stats > ${library}.oxoG_metrics_lib.log
  """
}

process wgs_metrics_lib {

  tag{ library }

  label 'picard'

  publishDir params.outdir, mode: 'copy', overwrite: true,
  saveAs: {filename ->
    if (filename.endsWith("log")) "pipeline/logs/$filename"
    else if (filename.endsWith("stats")) "qc/data/$filename" }
	
  input:
  set library, indiv, bamfile from remove_duplicates_wgs_metrics_lib_ch

  output:
  file "*.log" into wgs_metrics_lib_log_ch
  file "*.stats" into wgs_metrics_lib_stats_ch
  
  script:
  def java_mem = new String( task.memory.toString() ).replace( " GB", "" )
  """
  java -Xmx${java_mem}g -Djava.io.tmpdir=./ \
    -jar \${EBROOTPICARD}/picard.jar CollectWgsMetrics \
    R=${params.fasta} \
    I=${bamfile} \
    O=${library}.wgs_metrics_lib.stats > ${library}.wgs_metrics_lib.log
  #java -Xmx${java_mem}g -Djava.io.tmpdir=./ \
  #  -jar \${EBROOTPICARD}/picard.jar CollectWgsMetricsWithNonZeroCoverage  \
  #  R=${params.fasta} \
  #  I=${bamfile} \
  #  CHART=chart.pdf \
  #  O=${library}.wgs_metrics_nonzero_lib.stats > ${library}.wgs_metrics_nonzero_lib.log
  """
}

process collect_alignment_metrics_library {

  tag{ library }

  label 'picard'

  publishDir params.outdir, mode: 'copy', overwrite: true,
  saveAs: {filename ->
    if (filename.endsWith("log")) "pipeline/logs/$filename"
    else if (filename.endsWith("stats")) "qc/data/$filename" }
	
  input:
  set library, indiv, bamfile from remove_duplicates_collect_alignment_metrics_library_ch

  output:
  file "*.log" into collect_alignment_metrics_library_log_ch
  file "*.stats" into collect_alignment_metrics_library_stats_ch
  
  script:
  def java_mem = new String( task.memory.toString() ).replace( " GB", "" )
  """
  java -Xmx${java_mem}g -Djava.io.tmpdir=./ \
    -jar \${EBROOTPICARD}/picard.jar CollectAlignmentSummaryMetrics \
    R=${params.fasta} \
    I=${bamfile} \
    O=${library}.collect_alignment_metrics.stats > ${library}.collect_alignment_metrics.log
  """
}

process bam_to_vcf_library {

  tag{ library }

  publishDir params.outdir, mode: 'copy', overwrite: true,
  saveAs: {filename ->
    if (filename.endsWith("log")) "pipeline/logs/$filename"
    else if (filename.endsWith("stats")) "qc/data/$filename" }
	
  input:
  set library, indiv, bamfile from remove_duplicates_bam_to_vcf_library_ch

  output:
  file "*.log" into bam_to_vcf_library_log_ch
  file "*.stats" into bam_to_vcf_library_stats_ch
  
  script:
  """
  samtools mpileup -vf ${params.fasta} ${bamfile} > ${library}.vcf.gz 2> ${library}_bam_to_vcf_library.log
  bcftools stats ${library}.vcf.gz > ${library}.vcf.stats
  """
}


/* Merge libraries from the same individuals
 */
remove_duplicates_ch
    .groupTuple(by: 1)
    .map{ s -> tuple( s[0], s[1][0], s[2] ) } 
    .set{ indiv_merge_bams_ch }

process merge_bams {

  tag{ indiv }

  publishDir params.outdir, mode: 'copy', overwrite: true,
  saveAs: {filename ->
    if( filename.endsWith(".bam")) "results/alignments/$filename"
    else if (filename.endsWith(".bai")) "results/alignments/$filename"
    else if (filename.endsWith("flagstat" )) "qc/data/$filename"
    else if (filename.endsWith("idxstats" )) "qc/data/$filename" }
    
  input:
  set librarys, indiv, bamfiles from indiv_merge_bams_ch

  output:
  set indiv, "*.bam" into(
    merged_bams_pmdtools_metrics_ch,
    merged_bams_bam_to_vcf_ch,
    merged_bams_insert_size_ch,
    merged_bams_wgs_metrics_ch,
    merged_bams_collect_alignment_metrics_ch,
    merged_bams_randfa_ch )
  file "*.flagstat" into merge_bams_flagstat_ch
  file "*.idxstats" into merge_bams_idxstats_ch
  
  script:
  def bams_string = bamfiles.join(" ")
  //def indiv = indiv[0]
  """
  samtools merge ${indiv}.merged.bam ${bams_string}
  samtools sort -@ ${task.cpus} -o ${indiv}.bam ${indiv}.merged.bam
  samtools index ${indiv}.bam
  samtools flagstat ${indiv}.bam > ${indiv}.bam.flagstat
  samtools idxstats ${indiv}.bam > ${indiv}.bam.idxstats
  rm ${indiv}.merged.bam
  """
}

process pmdtools {

  tag{ indiv }

  errorStrategy 'ignore'

  publishDir params.outdir, mode: 'copy', overwrite: true,
  saveAs: {filename ->
    if (filename.endsWith("log")) "pipeline/logs/$filename"
    else if (filename.endsWith("stats")) "qc/data/$filename" }
	
  input:
  set indiv, bamfile from merged_bams_pmdtools_metrics_ch

  output:
  file "*.log" into pmd_metrics_log_ch
  file "*.stats" into pmd_metrics_stats_ch
  
  script:
  """
  samtools view -q30 ${bamfile} | pmdtools --first \
  2> ${indiv}.pmdtools.log
  cp ${indiv}.pmdtools.log ${indiv}.pmdtools.stats
  """
}

process wgs_metrics {

  tag{ indiv }

  label 'picard'

  publishDir params.outdir, mode: 'copy', overwrite: true,
  saveAs: {filename ->
    if (filename.endsWith("log")) "pipeline/logs/$filename"
    else if (filename.endsWith("stats")) "qc/data/$filename" }
	
  input:
  set indiv, bamfile from merged_bams_wgs_metrics_ch

  output:
  file "*.log" into wgs_metrics_log_ch
  file "*.stats" into wgs_metrics_stats_ch
  
  script:
  def java_mem = new String( task.memory.toString() ).replace( " GB", "" )
  """
  java -Xmx${java_mem}g -Djava.io.tmpdir=./ \
    -jar \${EBROOTPICARD}/picard.jar CollectWgsMetrics \
    R=${params.fasta} \
    I=${bamfile} \
    O=${indiv}.wgs_metrics.stats > ${indiv}.wgs_metrics.log
  #java -Xmx${java_mem}g	-Djava.io.tmpdir=./ \
  #  -jar \${EBROOTPICARD}/picard.jar CollectWgsMetricsWithNonZeroCoverage  \
  #  R=${params.fasta} \
  #  I=${bamfile} \
  #  CHART=chart.pdf \
  #  O=${indiv}.wgs_metrics_nonzero.stats > ${indiv}.wgs_metrics_nonzero.log
  """
}

process collect_alignment_metrics {

  tag{ indiv }

  label 'picard'

  publishDir params.outdir, mode: 'copy', overwrite: true,
  saveAs: {filename ->
    if (filename.endsWith("log")) "pipeline/logs/$filename"
    else if (filename.endsWith("stats")) "qc/data/$filename" }
	
  input:
  set indiv, bamfile from merged_bams_collect_alignment_metrics_ch

  output:
  file "*.log" into collect_alignment_metrics_log_ch
  file "*.stats" into collect_alignment_metrics_stats_ch
  
  script:
  def java_mem = new String( task.memory.toString() ).replace( " GB", "" )
  """
  java -Xmx${java_mem}g -Djava.io.tmpdir=./ \
    -jar \${EBROOTPICARD}/picard.jar CollectAlignmentSummaryMetrics \
    R=${params.fasta} \
    I=${bamfile} \
    O=${indiv}.collect_alignment_metrics.stats > ${indiv}.collect_alignment_metrics.log
  """
}

/* bamtovcf
   Convert bam file to vcf format via pileup
*/
process bam_to_vcf {

  tag{ indiv }

  publishDir params.outdir, mode: 'copy', overwrite: true,
  saveAs: {filename ->
    if( filename.endsWith(".vcf.gz") ) "results/vcf/$filename"
    else if (filename.endsWith("log")) "pipeline/logs/$filename"
    else if (filename.endsWith("stats")) "qc/data/$filename" }
	
  input:
  set indiv, bamfile from merged_bams_bam_to_vcf_ch

  output:
  set indiv, "*.vcf.gz" into bam_to_vcf_vcf_to_consensusfa_ch, bam_to_vcf_hets_ch
  file "*.log" into bam_to_vcf_log_ch
  file "*.stats" into bam_to_vcf_stats_ch
  
  script:
  """
  samtools mpileup -vf ${params.fasta} ${bamfile} > ${indiv}.vcf.gz 2> ${indiv}_bam_to_vcf.log
  bcftools stats ${indiv}.vcf.gz > ${indiv}.vcf.stats
  """
}

/* vcf_to_consensusfa
 */
process vcf_to_fasta {

  tag{ indiv }

  publishDir params.outdir, mode: 'copy', overwrite: true,
  saveAs: {filename ->
	    if( filename.endsWith(".fa")) "results/consensus/$filename"
	else if (filename.endsWith("log")) "pipeline/logs/$filename" }
	    
  input:
  set indiv, vcffile from bam_to_vcf_vcf_to_consensusfa_ch

  output:
  set indiv, "*.fa" into vcf_to_consensusfa_ch
  file "*.log" into vcf_to_consensusfa_log_ch

  script:
  """
  vcftools consensus \
   -I \
   -f \
   ${params.fasta} \
   ${vcffile} > ${indiv}.fa 2> ${indiv}_vcf_to_fasta.log
  """
}

/* het/diploid calls
 */
process randfa {

  tag{ indiv }

  publishDir params.outdir, mode: 'copy', overwrite: true,
  saveAs: {filename ->
    if( filename.endsWith(".randfa")) "results/rand/$filename"
	else if (filename.endsWith("log")) "pipeline/logs/$filename" }
		
  input:
  set indiv, bamfile from merged_bams_randfa_ch 

  output:
  set indiv, "*.randfa" into randfa_ch
  file "*.log" into randfa_log_ch

  script:
  """
  htsbox pileup \
  -R \
  -q 30 \
  -Q 30 \
  -l 35 \
  -s 1 \
  ${bamfile} > ${indiv}.randfa 2>${indiv}_randfa.log
  """
}

/* multiqc */
process multiqc {

  tag{ "multiqc" }

  publishDir "qc", mode: 'copy', overwrite: true

  input: file remove_duplicates_metrics_file from remove_duplicates_metrics_ch.collect()
  input: file remove_duplicates_flagstat_file from remove_duplicates_flagstat_ch.collect()
  input: file remove_duplicates_idxstats_file from remove_duplicates_idxstats_ch.collect()
  input: file collect_alignment_metrics_library_file from collect_alignment_metrics_library_stats_ch.collect()
  input: file collect_alignment_metrics_file from collect_alignment_metrics_stats_ch.collect()
  input: file merge_bams_flagstat_file from merge_bams_flagstat_ch.collect()
  input: file merge_bams_idxstats_file from merge_bams_idxstats_ch.collect()
  input: file vcf_stats_library_file from bam_to_vcf_library_stats_ch.collect()
  input: file vcf_stats_file from bam_to_vcf_stats_ch.collect()
  input: file wgs_metrics_stats_file from wgs_metrics_stats_ch.collect()
  input: file wgs_metrics_lib_stats_file from wgs_metrics_lib_stats_ch.collect()

  input: file oxoG_metrics_lib_stats_file from oxoG_metrics_lib_stats_ch.collect()
  
  script:
  """
  cd ${workflow.launchDir}
  multiqc -f -c ${params.multiqc_config} -o qc qc/data
  cd -
  """
}