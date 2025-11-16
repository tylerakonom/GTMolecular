# GTMolecular
*Post-publication analysis of patients with metastatic colorectal cancer pre/post treatment*

## First steps
*Identify how these data are presented and what they represent*

Data is available at NCBI BioProject accession number PRJNA714799 and represents cell-free DNA (cfDNA) isolated from the plasma of patients with metastatic colorectal cancer (CRC).

Sequencing was performed by targeted deep sequencing for two multiple gene panels (hybrid capture) that targeted genes frequently mutated in CRC. Paired-end reads were generated in fastq format.

## Create an analysis plan
*For this kind of analysis, I plan to preprocess reads the same way that researchers performed their preprocessing.*

### Data retrieval

1. Download raw reads of a single patient pre/post treatment (Fastq)

*I picked CTDC33 as my sample to work with. This patient had baseline, first response, and PBMC samples all sequenced with the same targeted gene panel*

2. Identify relevant samples/runs in the NCBI Run Selector to download. Patient CTDC33 has four relevant samples:

* CTDC33 (Baseline) (SRR13973737)
* CTDC33-2 (1st response) (SRR13973738)
* CTDC33-3 (PBMCs) (SRR13973947)

3. Downloaded to local machine using [SRA toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) by navigating to local directory and running:

	fastq-dump [Run number]

### Preprocessing
*I don't have access to enough personal computing power to perform this analysis, so I initially used Google Colab. I ran into an issue with CPU throughput during alignment, so I moved my processing to [CU Boulder's Research Computing](https://www.colorado.edu/rc/).

1. Deduplicate reads based off of random barcodes on P7 index sites
* The authors don't elaborate on this step outside of "used in-house scripts". They also don't give much information on what step the publically-available dataset was uploaded at. However, the manuscript says that the reads are paired-end, and there is only a single file per sample. For this reason, I expect that the files are uploaded after filtering and adapater removal.

2. Map (Burrows-Wheeler Aligner ["mem" algorithm] [v0.7.17]) to appropriate reference genome (GRCh38)
* Instead of building indexes from scratch using bwa, I downloaded the necessary index file from NCBI directly.
* -M tag required for GATK variant calling
*This step proved to be too much to run on Google Colab. The aligner that the authors used is CPU limited, and Google Colab only allows the use of a single CPU. Processing 3 samples would've taken ~72 hours, so I downloaded the fastq files to my local machine and uploaded them to RC at CU Boulder. I no longer have access to software installation, so I performed the next steps using the CURC-provided modules. Documentation for CURC can be found [here](https://curc.readthedocs.io/en/latest/index.html).*
* Adapted previously used scripts to submit alignment jobs to CURC using 24 cores on 1 thread. [Alignment script](https://github.com/tylerakonom/GTMolecular/blob/main/shell_scripts/alignReads.sh) was uploaded to scratch space, and jobs were queued using [this script](https://github.com/tylerakonom/GTMolecular/blob/main/shell_scripts/run_alignReads.sh) on the command line.

3. Convert to BAM files and sorted (Picard [v2.27.5])
* Run using a compile node and command line on the alpine cluster at CU Boulder. Loaded picard module using:

$module load picard

* Performed read sorting on the command line with:

$ java -jar $PICARD SortSam INPUT={filename}.sam OUTPUT={filename}_sorted.bam SORT_ORDER=coordinate

* Generated alignment metrics with samtools (v1.16.1) on the command line using:

$ samtools flagstat {filename}.sam > {filename}_alignment.txt

*Alignment metrics can be found [here](https://github.com/tylerakonom/GTMolecular/tree/main/alignment_metrics).*

* Perform routine deduplication using picard on the command line using:

$ java -jar $PICARD MarkDuplicates INPUT={filename}_sorted.bam OUTPUT={filename}_dedup.bam METRICS_FILE={filename}_dedup_metrics.txt

*Deduplication metrics can be found [here](https://github.com/tylerakonom/GTMolecular/tree/main/deduplication_metrics)*









### Data analysis

1. Perform variant calling

* Local realignment and base quality score recalibration (GATK [v4.1.0.0])
* Generate pileup files with SAMttools (mpileup)
* Variant calling (Varscan12)
* Annotation (SnpEff)
* Variant information added from ClinVar and COSMIC

2. Determine analysis parameters

* Filter out variants with less than 10% variant allele frequency (VAF) or less than 10 variant reads
    


