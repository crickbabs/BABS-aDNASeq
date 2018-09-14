
# BABS-aDNASeq

A Nextflow pipeline script for processing aDNASeq samples.

The pipeline was written by The Bioinformatics & Biostatistics Group in collaboration with The Ancient Genomics Lab @ The Francis Crick Institute.

nextflow: http://www.nextflow.io
nextflow-quickstart: http://www.nextflow.io/docs/latest/getstarted.html#get-started

### The Analysis Pipeline Flow

![alt unistrap](docs/BABS-aDNASeq_map.png) 

### Quick Start

To run an BABS-aDNASeq analysis you will need to complete the following steps. These are explained in more detail further down.

1) Obtain BABS-aDNASeq files from GitHub.
2) Install/load nextflow-0.32.0 or higher.
3) Configure reference genome file paths (genome.yml).
4) Configure environment profile if running software via a module system.
5) Create a sample design file.
6) Run nextflow pipeline.

### Get BABS-aDNASeq Files

To obtain BABS-aDNASeq files run the following git command.

   git clone https://github.com/crickbabs/BABS-aDNASeq

BABS-aDNASeq.nf		The Nextflow script.
BABS-aDNASeq			Wrapper script to run an analysis.
nextflow.config			Main BABS-aDNASeq config file.
conf/babs_profile.config	Profile configuration for running the script @ The Crick.
conf/genomes.config		Genomes configuration file for defining reference data.
conf/multiqc_config.yml		Multiqc configuration used to generate integrated QC report.


### Load Nextflow Module

If you are working within a module environment such as that at The Crick, load the nextflow module.

   module purge
   module load nextflow/0.30.2

### Sample Design File

Fastq files are specified in a csv design file with the following columns.

      column 1 : Individual ID
      column 2 : Sequencing library ID
      column 3 : full path to fastq file R1
      column 4 : full path to fastq file R2

### Running an aDNASeq-ByBABS Analysis

    BABS-aDNASeq --outdir ./ --design design.csv --profile babs --genome hg19 --resume

### Output Directories & Files

### Flow Details

#### Merge fastqs with the same library ID

#### Adapter trimming with SeqPrep

Adapter trimming and paired-end overlap consensus building. Only the overlap is saved here. Non-overlapping read-pairs are discarded.

https://github.com/jstjohn/SeqPrep

#### BWA

Consensus overlaps are aligned to the specified reference using BWA. BAM files with read groups are created.

https://github.com/lh3/bwa/

#### Duplicate Removal

Duplicate alignments are removed using Picard.

#### Variant Calling

VCFs are created using samtools mpileup. QC metrics are produced using bcftools stats.

#### Consensus Fasta

Ambiguity encoded consensus fasta files are produced using vcftools consensus.

#### Random Fasta

Random allele fasta files are produced using htsbox pileup -R.

https://github.com/lh3/htsbox

#### Merge BAM files

BAM files from the same individual are merged using samtools merge. Varient calling and QC ae carried out at both the library and individual level. 

#### QC

Alignment QC is assessed using pmdtools and CollectWgsMetrics, CollectWgsMetricsWithNonZeroCoverage, CollectOxoGMetrics & CollectAlignmentSummaryMetrics from Picard. A QC report is generated using multiqc.

     https://github.com/pontussk/PMDtools
     https://github.com/broadinstitute/picard
     https://github.com/ewels/MultiQC

#### Credits

The BABS-aDNASeq nextflow pipeline was written and developed by Philip East & Pontus Skoglund.

The Bioinformatics & Biostatistics Group (BABS) @ The Francis Crick Institute.
Ancient Genomics @ The Francis Crick Institute.

#### Licence

This project is licensed under the MIT License - see the LICENSE.md file for details.