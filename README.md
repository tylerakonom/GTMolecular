# GTMolecular
*Post-publication analysis of a patient ctDNA sample pre/post treatment for colorectal cancer as an exercise for GT Molecular*

## First steps
*Identify how these data are presented and what they represent*

* Data is available at NCBI BioProject accession number PRJNA714799 and represents cell-free DNA (cfDNA) isolated from the plasma of patients with metastatic colorectal cancer (CRC).
* Sequencing was performed by targeted deep sequencing for two multiple gene panels (hybrid capture) that targeted genes frequently mutated in CRC. Paired-end reads were generated in fastq format.

### Create an analysis plan
*To start, I plan to process and analyze reads the same way that the original authors did. This will give me an in-depth idea of how to improve on their methods for future analysis plans.*

The steps I plan to take to generate visualized results are as follows:
1. Retrieve raw reads
	* Samples will be received in the form of fastq files, representing millions of raw reads from the sequencer. Without processing, there is no way to contextualize these data.
2. Preprocess samples
	* Remove adapter sequences: adapter sequences are attached to each DNA fragment during NGS library preparation. These sequencers are not representative of the host's genome, and can cause issues when attempting to align reads to the genome.
	* Aligning: cDNA reads are compared to a reference genome, then assigned a place in the genome that they likely came from.
3. Base quality score recalibration
	* During NGS, bases are assigned a quality score based off of the probability that the base call is correct and are regulary higher than they should be.
	* Using machine learning, a reference genome, and a database of known variants we recalibrate these quality scores more accurately assess quality scores.
4. Perform variant calling
	* GATK and a reference genome are used to identify variations in our mapped reads compared to a refernce genome.
	* Variations are filtered out by setting accepted QC thresholds, generating a final vcf with variants of interest.

## Data processing

### Data retrieval
1. Download raw reads of a single patient pre/post treatment (Fastq)
	* I picked CTDC33 as my sample to work with. This patient had baseline, first response, and PBMC samples all sequenced with the same targeted gene panel

2. Identify relevant samples/runs in the NCBI Run Selector to download. Patient CTDC33 has four relevant samples:
	* CTDC33 (Baseline) (SRR13973737)
	* CTDC33-2 (1st response) (SRR13973738)
	* CTDC33-3 (PBMCs) (SRR13973947)

3. Created a conda environment in Google Colab and installed [SRA toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit), then downloaded files using:
	```bash
	$ fastq-dump [Run number]
	```

### Preprocessing
*I don't have access to enough personal computing power to perform this analysis, so I initially used Google Colab and [this notebook](https://github.com/tylerakonom/GTMolecular/blob/main/colab_notebook/GTMolecular.ipynb). I ran into an issue with CPU throughput during alignment, so I moved my processing to [CU Boulder's Research Computing](https://www.colorado.edu/rc/).*

1. Deduplicate reads based off of random barcodes on P7 index sites
	* The authors don't elaborate on this step outside of "used in-house scripts", but I performed deduplication as part of a post-alignment step. They also don't give much information on what step the publically-available dataset was uploaded at. However, the manuscript says that the reads are paired-end, and there is only a single file per sample. For this reason, I expect that the files are uploaded after filtering and adapater removal.

2. Map using Burrows-Wheeler Aligner ("mem" algorithm) (v0.7.17) to appropriate reference genome (GRCh38).
	* Instead of building indexes from scratch using bwa, I downloaded the necessary index file from NCBI directly using:
		```bash
		$ wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_analysis_set.fna.bwa_index.tar.gz
		$ gunzip -q GCA_000001405.15_GRCh38_full_analysis_set.fna.bwa_index.tar.gz
		$ tar -xvf GCA_000001405.15_GRCh38_full_analysis_set.fna.bwa_index.tar
		```
	* -M tag required for GATK variant calling
	*This step proved to be too much to run on Google Colab. The aligner that the authors used is CPU limited and Google Colab only allows the use of a single CPU. Processing 3 samples would've taken ~72 hours, so I downloaded the fastq files to my local machine and uploaded them to RC at CU Boulder. I no longer have access to software installation, so I performed the next steps using the CURC-provided modules. Documentation for CURC can be found [here](https://curc.readthedocs.io/en/latest/index.html).*
	* Adapted previously used scripts to submit alignment jobs to CURC using 24 cores on 1 thread. [Alignment script](https://github.com/tylerakonom/GTMolecular/blob/main/shell_scripts/alignReads.sh) was uploaded to scratch space, and jobs were queued using [this script](https://github.com/tylerakonom/GTMolecular/blob/main/shell_scripts/run_alignReads.sh) on the command line.

3. Convert to BAM files and sorted (Picard) (v2.27.5)
	* Run using a compile node and command line on the alpine cluster at CU Boulder. Loaded picard module using:
		```bash
		$ module load picard
		```
	* Performed read sorting on the command line with:
		```bash
		$ java -jar $PICARD SortSam INPUT={filename}.sam OUTPUT={filename}_sorted.bam SORT_ORDER=coordinate
		```
	* Generated alignment metrics with samtools (v1.16.1) on the command line using:
		```bash
		$ samtools flagstat {filename}.sam > {filename}_alignment.txt
		```
		*Alignment metrics can be found [here](https://github.com/tylerakonom/GTMolecular/tree/main/alignment_metrics).*

	* Perform routine deduplication using picard on the command line using:
		```bash
		$ java -jar $PICARD MarkDuplicates INPUT={filename}_sorted.bam OUTPUT={filename}_dedup.bam METRICS_FILE={filename}_dedup_metrics.txt
		```
		*Deduplication metrics can be found [here](https://github.com/tylerakonom/GTMolecular/tree/main/deduplication_metrics)*


### Variant calling

1. Local realignment and base quality score recalibration (BQSR) (GATK) (v4.1.0.0)
	* Generated reference dictionary, fasta index, and bam index for GATK:
		```bash
		$ java -jar picard CreateSequenceDictionary R=GCA_000001405.15_GRCh38_full_analysis_set.fna O=GCA_000001405.15_GRCh38_full_analysis_set.fna.dict
		$ samtools faidx GCA_000001405.15_GRCh38_full_analysis_set.fna
		$ samtools index {filename}_dedup.bam
		```
	* Downloaded "known variant" file and index following GATK best practices located on the [Broad Institute's GitHub](https://github.com/gatk-workflows/gatk4-data-processing/blob/master/processing-for-variant-discovery-gatk4.wdl), with the file bucket located [here](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?pli=1&prefix=&forceOnObjectsSortingFiltering=false).
	* Performed BQSR using [this](https://github.com/tylerakonom/GTMolecular/blob/main/shell_scripts/bqsr.sh) script queued in the command line with [this](https://github.com/tylerakonom/GTMolecular/blob/main/shell_scripts/run_bqsr.sh) script. GATK was called twice, the first time to generate a BQSR table, and the second time was to apply the new scores to .bam files.

2. Perform variant calling
	* Variant calling (GATK) performed using [this script](https://github.com/tylerakonom/GTMolecular/blob/main/shell_scripts/caller.sh) and jobs were queued using [this script](https://github.com/tylerakonom/GTMolecular/blob/main/shell_scripts/run_caller.sh). GATK was called 5 total times. (1) I generated a .vcf of "raw variants", then (2,3) I selected SNPs and INDELs from this file and separated them into their own files. Finally (4,5) I filtered the variants for various QC metrics to remove false positives.
	* Annotation (SnpEff) to the NCBI GRCh38 refseq database was performed on the command line using:
		```bash
		$ snpEff download -v GRCh38.mane.1.2.refseq
		$ snpEff -v GRCh38.mane.1.2.refseq {filename}_filtered_snps.vcf > {filename}_filtered_snps.ann.vcf
		```
	* All vcf files generated are located [here](https://github.com/tylerakonom/GTMolecular/blob/main/variants/).
	* Data were visualized using IGV.


## Results

### Pipeline troubleshooting

Since I only have limited access to previous research computing resources, I had to adapt my approach from my previous experience with CURC to Google Colab. Once I was forced to pivot back to CURC due to throughput issues I was forced to use the tools that are provided instead of the exact matches to the tools the author's used in the manuscript. This means that, even though I arrived at visualized results, I expect them to deviate from the original author's. Another consequence of having to pivot multiple times is that I ran out of time to have the final step of sample 3 run completely. I have a finalized .bam file for sample CTDC33_3, but was unable to generate a final SNP or INDEL file.

### Sample discrepancies

The samples uploaded to NCBI's Bioproject repository were not named in a consistent way, and it seems that the authors had a "raw" sample ID that differed greatly from the sample IDs that they used in their bioinformatics pipeline (supplementary table 2). I ended up downloading samples named CTC333, CTC333-2, and PBMC_CTC333 expecting that they were samples CTDC33 (baseline), CTDC33-2 (1st response), and CTDC33-3 (PBMCs) analyzed with "Version 3" of their NGS panel. After processing these data and visualizing them in IGV I found that NCBI sample CTC333 had very little coverage (confirmed by alignment metrics and [visualized in IGV](https://github.com/tylerakonom/GTMolecular/blob/main/igv_snapshots/1v2IGVSnapshot.png)), leading me to the conclusion that it was either a sample the authors weren't able to include in their final manuscript or something else entirely. Furthermore, the alignment metrics for each of the three samples don't match the supplementary table provided in the manuscript.

Unfortunately, these challenges result in a highly limited context in which we can draw conclusions. We can't determine whether or not this sample was from ctDNA and representative of colorectal cancer cells, or if this sample was from the PBMCs and representative of the patient's genome.

### Results from CTDC33_2 for *KRAS*

In sample CTDC33_2 (NCBI Bioproject sample CTC333-2), [two SNPs were identified](https://github.com/tylerakonom/GTMolecular/blob/main/igv_snapshots/kras.png) in protein coding (exon) sequences of the gene *KRAS* through variant calling analysis.

Depending on the context with which we wanted to present these results, I can imagine a couple of different ways to market them:

	* Advertising graphic for NGS customized pipeline services/library prep kit: I invision a simplified diagram of DNA with two obvious red regions with a tag line like "Through the use of GT Molecular's customized bioinformatics pipeline/GT Molecular's NGS library prep kit we were able to perform a deep sequencing of 16 genes sensitive enough to detect single nucleotide polymorphisms in *KRAS*, a gene highly relevant to the study of colorectal cancer" underneath it. Another option would be an annotated and worked-up diagram (similar to the [IGV snapshot](https://github.com/tylerakonom/GTMolecular/blob/main/igv_snapshots/kras.png) in a more digestible format) with various comparative metrics in a pared-down table underneath.

	* Comparison between this analysis and a GT Molecular customized bioinformatics pipeline: This analysis was performed using a very basic analysis method and outdated tools to match the author's pipeline as closely as possible within my time constraints. In my experience, academic researchers don't have the luxury to pursue the most cutting edge analysis methods due to a lack of time, budget, or expertise. This analysis could be performed again using an updated pipeline using the most relevant and powerful tools, expecting an increase in sensitivity. These results would then be contrasted to highlight the benefits of using a GT Molecular service.

## Future directions

	* Step-wise improvements: This process was performed using a basic variant calling workflow, and has tons of room to be improved. There are newer tools, more specific or relevant databases, and settings/thresholds that would require more research and time than I have to improve on this analysis drastically.
	* Snakemake/Nextflow: The first time an anlaysis is performed is always the hardest, and follow-up workflows can be heavily simplified or automated entirely. The most comprehensive solution I can think of to remedy this is to create a Snakemake or Nextflow script with a basic GUI and user manual. This would make this analysis much more accessible and theoretically marketable.


Thank you so much for your time reading through this analysis, and I look forward to discussing any feedback or critique you have for my processes!
