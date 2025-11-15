# GTMolecular
*Post-publication analysis of patients with metastatic colorectal cancer pre/post treatment*

## First steps
*Identify how these data are presented and what they represent*

Data is available at NCBI BioProject accession number PRJNA714799 and represents cell-free DNA (cfDNA) isolated from the plasma of patients with metastatic colorectal cancer (CRC).

Sequencing was performed by targeted deep sequencing for two multiple gene panels (hybrid capture) that targeted genes frequently mutated in CRC. Paired-end reads were generated in fastq format.

## Create an analysis plan
*For this kind of analysis, I plan to preprocess reads the same way that researchers performed their preprocessing.*

### Data preprocessing

1. Download raw reads of a single patient pre/post treatment (Fastq)

*I picked CTDC33 as my sample to work with. This patient had baseline, first response, disease progression at death, and PBMC samples all sequenced with the same targeted gene panel*

Identify relevant samples/runs in the NCBI Run Selector to download

Patient CTDC33 has four relevant samples:

* CTDC33 (Baseline) (SRR13973737)
* CTDC33-2 (1st response) (SRR13973738)
* CDTC33-3 (PBMCs) (SRR13973947)

Downloaded to local machine using [SRA toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) by navigating to local directory and running:

	fastq-dump [Run number]

2. Preprocess/prepare data
    
* Map (Burrows-Wheeler Aligner [v0.7.10] ["mem" algorithm]) to appropriate reference genome (hg19)
* Convert to BAM files and sorted (SAMtools [v1.1])


### Data analysis

3. Perform variant calling

* Local realignment and base quality score recalibration (GATK [v4.1.0.0])
* Generate pileup files with SAMttools (mpileup)
* Variant calling (Varscan12)
* Annotation (SnpEff)
* Variant information added from ClinVar and COSMIC

4. Analysis parameters

* Filter out variants with less than 10% variant allele frequency (VAF) or less than 10 variant reads
    


