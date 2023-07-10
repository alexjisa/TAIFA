> **TAIFA is a workflow for conducting Whole-Genome Sequencing (WGS) of bacteria**

# Introduction
>
WGS provides a global view of the entire genome. The flexible, scalable and fast nature of Next-Generation Sequencing (NGS) techniques as well as the dropping costs in recent years opens the door to a powerful tool for genomics research. Sequencing, de novo assembly and annotation of bacterial genomes has great potential for clinical, industrial or basic science applications.
>
# Context
>
Bacterial genomic DNA is obtained using an extraction and purification kit according to the instructions of the commercial company.
>
Libraries are prepared using Illumina DNA Prep technology (previously known as Nextera DNA Flex).
>
Quantitation and quality control of libraries is performed using the PicoGreen assay and the High Sensitivity DNA assay on the 2100 Bioanalyzer System, respectively.
>
Results are from Illumina iSeq100 paired-end sequencing run, 2 × 150 bp. Around 17 hours with a maximum output of 1.2 Gb.
>
R1.fastq.gz, R2.fastq.gz, Undetermined_R1.fastq.gz and Undetermined_R2.fastq.gz are automatically generated from the base-calling (bci).

# Workflow summary
1. Sequencing quality control and trimming
2. De novo genome assembly and evaluation 
3. Mapping 
4. Annotation
5. Pathogenicity prediction 
  
# Workflow
## 1. Sequencing quality control and trimming
*1.1. First FastQC*
>
Description:
>
Basic quality control checks of the raw reads by browsing an html file containing the main quality parameters.  
>
Link: https://github.com/s-andrews/FastQC/blob/master/README.md 
```{bash}
/path/to/fastqc path/to/*.fastq.gz -o /path/to/fastqc_output
```
*1.2. Trimmomatic*
>
Description:

>
Link: https://github.com/usadellab/Trimmomatic/blob/main/README.md
```{bash}
cd /path/to/
java -jar /path/to/trimmomatic-0.39.jar PE -phred33 R1.fastq.gz R2.fastq.gz trim_R1.fastq.gz Undetermined_R1.fastq.gz trim_R2.fastq.gz Undetermined_R2.fastq.gz CROP:150 HEADCROP:1 AVGQUAL:20 SLIDINGWINDOW:4:20
```
*1.3. Second FastQC*
>
Description:
>
This quality control allows us to check if we are satisfied with the cleaning of the reads before proceeding to the next steps.
```{bash}
/path/to/fastqc path/to/trim_* -o /path/to/fastqc_output
```
## 2. De novo genome assembly and evaluation 
*2.1. A5-miseq*
>
Description:
Pipeline for assembling DNA sequence data generated on the Illumina sequencing platform (homozygous haploid genomes and reads greater than 80 pb).
>
Easy, fast and the best results.
>
Link: https://sourceforge.net/p/ngopt/wiki/A5PipelineREADME/ ; DOI:10.1093/bioinformatics/btu661
```{bash}
/path/to/a5_pipeline.pl --threads=4 R1.fastq R2.fastq a5_output 
```

*2.2. QUAST*
>
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
>
>
Description:
Evaluation according to assembly continuity (number of bp, number of contigs or N50) and other information (like GC%, largest contig or mismatches).
>
Link: https://quast.sourceforge.net/docs/manual.html
```{bash}
conda activate quast_env
quast.py a5_output.contigs.fasta -o quast_a5_output -l a5
````

*2.3. BUSCO*
>
Description:
Evaluation according to assembly content.
>
Link: https://busco.ezlab.org/v3/
```{bash}
module load Anaconda3/4.4.0
module load augustus
run_BUSCO.py -c 4 -i a5_output.contigs.fasta -l /path/to/bacteria_odb9 -o busco_a5 -m geno 
```

## 3. Mapping
*3.1. BWA, SAMtools*
>
Description:
>
Link:
```{bash}
bwa index mygenome.contigs.fasta
bwa mem mygenome.contigs.fasta trim_R1_concatenado.fastq.gz trim_R2_concatenado.fastq.gz > result.sam
samtools view -b -o result.bam result.sam
Tras ejecutar este .sh procedemos a comprobar si todo ha ido bien en result.bam:
samtools view -H result.bam
samtools view result.bam | tail
rm -fr result.sam
samtools sort result.bam -o result_sort.bam
samtools index result_sort.bam
samtools faidx result_sort.bam
```

*3.2. IGV*
>
Description:
>
Link: https://software.broadinstitute.org/software/igv/

*3.3. Qualimap*
>
Description:
>

>
Link: http://qualimap.conesalab.org/

## 4. Annotation
*4.1. Structural annotation: Prokka*
>
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
>
>
Description:
>
It combines the use of Prodigal (location of coding genes) and Infernal (location of non-protein coding genes) to generate as output several files (faa, ffn, gff or gbk, among others) that will be required for functional annotation.
>
Link: https://github.com/tseemann/prokka/blob/master/README.md
```{bash}
conda activate prokka_env
prokka a5_output.contigs.fasta --cpus 6 --outdir prokka_output
```

*4.2. Functional annotation: EggNOG-mapper*
>
Description:
>
EggNOG-mapper uses DIAMOND and HMM to functionally annotate using groups of orthologues at different taxonomic levels and phylogenies from the eggNOG database.
>
Its public online resource (http://eggnog-mapper.embl.de/) allows uploading up to 100,000 CDS in FASTA format, more than enough to annotate a complete bacterial genome.
>
The annotation is one of the most complete and the execution time is very affordable.
>
Link: https://github.com/eggnogdb/eggnog-mapper

## 5. Pathogenicity prediction
>
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
>
>
Description:
>
Preliminary approaches to determine the possible pathogenicity of bacteria are necessary in any kind of research. The characterisation of virulence factors (VFs), proteins enable to cause infection, from WGS data allows you to get an idea of the microorganism you are working with.
>
Link: https://github.com/bbuchfink/diamond/blob/master/README.md ; http://www.mgc.ac.cn/cgi-bin/VFs/v5/main.cgi
```{bash
conda activate diamond_env
diamond blastp --db /home/alejandro_jimenez/vfdb -q a5_output_contigs.fasta -o vfdb_results.tsv --query-cover 97 --id 50
```

# Credits
>
Done by Alejandro Jimenez-Sanchez for his MSc dissertation under the supervision of Dr. Alvaro Polonio.
>
The bioinformatics tools included in this repository were chosen based on the execution time and the quality of the results generated by each of them according to the bacterial genomes we worked with. Different programs were evaluated, mainly assemblers and annotators, which are not listed. It is recommended to test alternatives to those proposed here as each project is unique. 
