> **TAIFA is a workflow for conducting Whole-Genome Sequencing (WGS) of bacteria**

# Introduction
>
WGS provides a global view of the entire genome. The flexible, scalable and fast nature of Next-Generation Sequencing (NGS) techniques as well as the dropping costs in recent years opens the door to a powerful tool for genomics research. Sequencing, de novo assembly and annotation of bacterial genomes has great potential for clinical, industrial or basic science applications.
>
# Context
>
Bacterial genomic DNA is obtained using an extraction and purification kit according to the instructions of the commercial company.
>
Libraries are prepared using Illumina DNA Prep technology (previously known as Nextera DNA Flex). Quantitation and quality control are performed before sequencing.
>
Results are from Illumina iSeq100 paired-end (PE) sequencing run, 2 × 150 bp. Around 17 hours with a maximum output of 1.2 Gb.
>
R1.fastq.gz, R2.fastq.gz, Undetermined_R1.fastq.gz and Undetermined_R2.fastq.gz are automatically generated from the base calls in the BCL files.

# Workflow summary
1. Sequencing quality control and trimming
2. De novo genome assembly and evaluation 
3. Mapping (optional)
4. Annotation
5. Pathogenicity prediction (optional)
  
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
Trimmomatic is a tool for trimming Illumina FASTQ data and removing adapters. In this case, the adapters are automatically removed in advance, so in the command what is done is: remove the last and first base, remove the read if the average quality is less than 20, and examines the read 4 by 4 bases and eliminates it if the average quality is lower than the given value.
>
Link: https://github.com/usadellab/Trimmomatic/blob/main/README.md
```{bash}
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
>
Pipeline for assembling DNA sequence data generated on the Illumina sequencing platform (homozygous haploid genomes and reads greater than 80 pb).
>
Easy, fast and the best results.
>
Link: https://sourceforge.net/p/ngopt/wiki/A5PipelineREADME/ ; DOI:10.1093/bioinformatics/btu661
```{bash}
/path/to/a5_pipeline.pl --threads=4 trim_R1.fastq trim_R2.fastq a5_output 
```

*2.2. QUAST*
>
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
>
>
Description:
>
Evaluation according to assembly continuity. QUAST provides various metrics that quantitatively describe the assembly quality, such as the number and size of contigs, N50 length (a measure of contig/scaffold length), genome fraction coverage, misassemblies, indels, and duplicated regions. These metrics help researchers understand the completeness and accuracy of the assembly.
>
Link: https://quast.sourceforge.net/docs/manual.html
```{bash}
conda activate quast_env
quast.py a5_output.contigs.fasta -o quast_a5_output -l a5
````

*2.3. BUSCO*
>
Description:
>
Evaluation according to assembly content. The purpose of BUSCO is to evaluate the completeness and accuracy of an assembly by comparing it to a set of highly conserved genes that are expected to be present in most species. In this case it utilizes a database of orthologous genes that are considered to be universally present as single-copy genes across bacteria. These genes are highly conserved, meaning they have minimal sequence variation and are maintained in a single copy within a genome. The output provides quantitative measures, such as the number and percentage of complete, fragmented, and missing orthologs.
>
Link: https://busco.ezlab.org/v3/
```{bash}
module load Anaconda3/4.4.0
module load augustus
run_BUSCO.py -c 4 -i a5_output.contigs.fasta -l /path/to/bacteria_odb9 -o busco_a5 -m geno 
```

## 3. Mapping (optional)
*3.1. BWA and SAMtools*
>
Description:
>
Construction of the index for the sequenced genome and mapping reads against this genome thanks to BWA index and MEM, respectively. Conversion of SAM file to BAM, ordering of BAM file by genomic position and index construction by samtools.
>
Link: https://github.com/lh3/bwa/blob/master/README.md ; https://github.com/samtools/samtools
```{bash}
bwa index a5_output.contigs.fasta
bwa mem a5_output.contigs.fasta trim_R1.fastq.gz trim_R2.fastq.gz > align.sam

samtools view -b -o align.bam align.sam
samtools sort align.bam -o align_sort.bam
samtools index align_sort.bam
```

*3.2. IGV*
>
Description:
>
It is an alignment viewer in which we load the resulting BAM file as well as the GFF file of the structural annotation which will be obtained later.
>
Link: https://software.broadinstitute.org/software/igv/

*3.3. Qualimap*
>
Description:
>
The BAM QC analysis mode evaluates the quality of the alignment data contained in the BAM file. It returns information on reads, coverage, alignment quality...
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
The annotation is one of the most complete (COG category,  GOs, EC, KEGG, CAZY or PFAMs) and the execution time is very affordable.
>
Link: https://github.com/eggnogdb/eggnog-mapper

## 5. Pathogenicity prediction (optional)
>
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
>
>
Description:
>
Preliminary approaches to determine the possible pathogenicity of bacteria are necessary in any kind of research. The characterisation of virulence factors (VFs), proteins enable to cause infection, from WGS data allows you to know you are doing your project safely and securely.

>
For this purpose, DIAMOND, a fast and sensitive protein aligner, is used to achieve ultra-fast alignments against the Virulence Factors DataBase (VFDB). The chosen parameters ensure a high sequence alignment identification.
>
Link: https://github.com/bbuchfink/diamond/blob/master/README.md ; http://www.mgc.ac.cn/cgi-bin/VFs/v5/main.cgi
```{bash
conda activate diamond_env
diamond blastp --db /home/alejandro_jimenez/vfdb -q a5_output_contigs.fasta -o vfdb_results.tsv --query-cover 97 --id 50
```
## Automated pipeline implementation

# 1) Control de calidad

## 1.1) Primer FastQC

```{bash}
# Ruta de la carpeta principal
origen="/home/alejandro_jimenez/seq"
# Recorrer todas las carpetas dentro de la carpeta principal
for carpeta in "$origen"/*; do
if [ -d "$carpeta" ]; then
  # Entrar en la carpeta "Alignment_1"
  cd "$carpeta/Alignment_1" || continue
  # Verificar si la carpeta "Fastq" existe
  if [ -d "Fastq" ]; then
    # Entrar en la carpeta "Fastq"
    cd "$carpeta/Alignment_1/Fastq" || continue
    # Ejecutar Fastqc
    mkdir fastqc
    /home/alejandro_jimenez/fastqc/fastqc *.fastq.gz -o "$carpeta"/Alignment_1/Fastq/fastqc
    # Mostrar mensaje de confirmacion
    echo "FastQC de las reads de la carpeta "$carpeta" realizado"
  fi
fi
done
```

## 1.2) Checkpoint: Trimmomatic

```{bash}
# Ruta de la carpeta principal
origen="/home/alejandro_jimenez/seq"

# Recorrer todas las carpetas dentro de la carpeta principal
for carpeta in "$origen"/*; do
if [ -d "$carpeta" ]; then
  # Entrar en la carpeta "Alignment_1"
  cd "$carpeta/Alignment_1" || continue
  # Verificar si la carpeta "Fastq" existe
  if [ -d "Fastq" ]; then
    # Entrar en la carpeta "Fastq"
    cd "$carpeta/Alignment_1/Fastq" || continue
    if [ -d "fastqc" ]; then
      echo "Ingresa el comando de Trimmomatic a ejecutar para las reads en "$carpeta":"
      
      #java -jar /data/home/alejandro_jimenez/trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 MG_S1_L001_R1_001.fastq.gz MG_S1_L001_R2_001.fastq.gz /home/alejandro_jimenez/seq/mg/Alignment_1/Fastq/fastqc/trim_MG_S1_L001_R1_001.fastq.gz Undetermined_S0_L001_R1_001.fastq.gz /home/alejandro_jimenez/seq/mg/Alignment_1/Fastq/fastqc/trim_MG_S1_L001_R2_001.fastq.gz Undetermined_S0_L001_R2_001.fastq.gz CROP:150 HEADCROP:1 AVGQUAL:20 SLIDINGWINDOW:4:20

      #java -jar /data/home/alejandro_jimenez/trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 MUS4_S1_L001_R1_001.fastq.gz MUS4_S1_L001_R2_001.fastq.gz /home/alejandro_jimenez/seq/mus4/Alignment_1/Fastq/fastqc/trim_MUS4_S1_L001_R1_001.fastq.gz Undetermined_S0_L001_R1_001.fastq.gz /home/alejandro_jimenez/seq/mus4/Alignment_1/Fastq/fastqc/trim_MUS4_S1_L001_R2_001.fastq.gz Undetermined_S0_L001_R2_001.fastq.gz CROP:150 HEADCROP:1 AVGQUAL:20 SLIDINGWINDOW:4:20

      #java -jar /data/home/alejandro_jimenez/trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 MUS7_S1_L001_R1_001.fastq.gz MUS7_S1_L001_R2_001.fastq.gz /home/alejandro_jimenez/seq/mus7/Alignment_1/Fastq/fastqc/trim_MUS7_S1_L001_R1_001.fastq.gz Undetermined_S0_L001_R1_001.fastq.gz /home/alejandro_jimenez/seq/mus7/Alignment_1/Fastq/fastqc/trim_MUS7_S1_L001_R2_001.fastq.gz Undetermined_S0_L001_R2_001.fastq.gz CROP:150 HEADCROP:1 AVGQUAL:20 SLIDINGWINDOW:4:20
      
      #java -jar /data/home/alejandro_jimenez/trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 MR_S1_L001_R1_001.fastq.gz MR_S1_L001_R2_001.fastq.gz /home/alejandro_jimenez/seq/mr/Alignment_1/Fastq/fastqc/trim_MR_S1_L001_R1_001.fastq.gz Undetermined_S0_L001_R1_001.fastq.gz /home/alejandro_jimenez/seq/mr/Alignment_1/Fastq/fastqc/trim_MR_S1_L001_R2_001.fastq.gz Undetermined_S0_L001_R2_001.fastq.gz CROP:150 HEADCROP:1 AVGQUAL:20 SLIDINGWINDOW:4:20

      read comando
      eval "$comando"
      # Mostrar mensaje de confirmacion
      echo "Trimmomatic de las reads de la carpeta "$carpeta" realizado"
    fi
  fi
fi
done
```

## 1.3) Segundo FastQC

```{bash}
# Ruta de la carpeta principal
origen="/home/alejandro_jimenez/seq"

# Recorrer todas las carpetas dentro de la carpeta principal
for carpeta in "$origen"/*; do
if [ -d "$carpeta" ]; then
  # Entrar en la carpeta "Alignment_1"
  cd "$carpeta/Alignment_1" || continue
  # Verificar si la carpeta "Fastq" existe
  if [ -d "Fastq" ]; then
    # Entrar en la carpeta "Fastq"
    cd "$carpeta/Alignment_1/Fastq" || continue
    if [ -d "fastqc" ]; then
      cd "$carpeta/Alignment_1/Fastq/fastqc" || continue
      # Ejecutar Fastqc
      /home/alejandro_jimenez/fastqc/fastqc trim_*
      # Mostrar mensaje de confirmacion
      echo "FastQC de las reads trimmeadas de la carpeta "$carpeta" realizado"
    fi
  fi
fi
done
```

# 2) Ensamblaje de novo y evaluación

## 2.1) A5-miseq

```{bash}
# Ruta de la carpeta principal
origen="/home/alejandro_jimenez/seq"

# Expresión regular para los archivos R1
patron_r1="trim_.*_S1_L001_R1_001.fastq"

# Expresión regular para los archivos R2
patron_r2="trim_.*_S1_L001_R2_001.fastq"

# Recorrer todas las carpetas dentro de la carpeta principal
for carpeta in "$origen"/*; do
if [ -d "$carpeta" ]; then
  # Entrar en la carpeta "Alignment_1"
  cd "$carpeta/Alignment_1" || continue
  # Verificar si la carpeta "Fastq" existe
  if [ -d "Fastq" ]; then
    # Entrar en la carpeta "Fastq"
    cd "$carpeta/Alignment_1/Fastq" || continue
    if [ -d "fastqc" ]; then
      cd "$carpeta/Alignment_1/Fastq/fastqc" || continue
      
      # Buscar archivos R1 que coincidan con el patrón en la carpeta principal y subcarpetas
      archivos_r1=$(find "$carpeta/Alignment_1/Fastq/fastqc" -regex ".*/$patron_r1")

      # Buscar archivos R2 que coincidan con el patrón en la carpeta principal y subcarpetas
      archivos_r2=$(find "$carpeta/Alignment_1/Fastq/fastqc" -regex ".*/$patron_r2")
      
      # Verificar si se encontraron archivos R1 y R2
      if [ -n "$archivos_r1" ] && [ -n "$archivos_r2" ]; then
      # Ejecutar el comando sobre cada par de archivos R1 y R2
        for archivo_r1 in $archivos_r1; do
          archivo_r2=$(echo "$archivo_r1" | sed 's/R1/R2/')
          nohup time -p /data/home/alejandro_jimenez/a5/bin/a5_pipeline.pl --threads=6 "$archivo_r1" "$archivo_r2" a5_output &
        done
        # Mostrar mensaje de confirmacion
        echo "Ensamblaje con A5 de las reads de la carpeta "$carpeta" realizado"
      else
      echo "No se encontraron archivos que coincidan con los patrones"
      fi
    fi
  fi
fi
done
```

## 2.2) QUAST

```{bash}
conda activate quast_env

# Ruta de la carpeta principal
origen="/home/alejandro_jimenez/seq"

# Recorrer todas las carpetas dentro de la carpeta principal
for carpeta in "$origen"/*; do
if [ -d "$carpeta" ]; then
  # Entrar en la carpeta "Alignment_1"
  cd "$carpeta/Alignment_1" || continue
  # Verificar si la carpeta "Fastq" existe
  if [ -d "Fastq" ]; then
    # Entrar en la carpeta "Fastq"
    cd "$carpeta/Alignment_1/Fastq" || continue
    if [ -d "fastqc" ]; then
      cd "$carpeta/Alignment_1/Fastq/fastqc" || continue
      quast.py a5_output.contigs.fasta -o quast_a5_output -l a5
    fi 
  fi
fi
done

conda deactivate
```

## 2.3) BUSCO

```{bash}
module load Anaconda3/4.4.0
module load augustus
run_BUSCO.py -c 4 -i a5_output.contigs.fasta -l /home/alejandro_jimenez/bacteria_odb9 -o busco_a5 -m geno 
```

# 3. Anotación estructural y funcional

## 3.1) Prokka

```{bash}
conda activate prokka_env

# Ruta de la carpeta principal
origen="/home/alejandro_jimenez/seq"

# Recorrer todas las carpetas dentro de la carpeta principal
for carpeta in "$origen"/*; do
if [ -d "$carpeta" ]; then
  # Entrar en la carpeta "Alignment_1"
  cd "$carpeta/Alignment_1" || continue
  # Verificar si la carpeta "Fastq" existe
  if [ -d "Fastq" ]; then
    # Entrar en la carpeta "Fastq"
    cd "$carpeta/Alignment_1/Fastq" || continue
    if [ -d "fastqc" ]; then
      cd "$carpeta/Alignment_1/Fastq/fastqc" || continue
      nohup time -p prokka a5_output.contigs.fasta --cpus 6 --outdir prokka_output &
    fi 
  fi
fi
done

cd /home/alejandro_jimenez/seq/mr/Alignment_1/Fastq/fastqc/spades_output
nohup time -p prokka contigs.fasta --cpus 6 --outdir prokka_output &

conda deactivate
```

## 3.2) eggNOGmapper

<http://eggnog-mapper.embl.de/> <https://github.com/eggnogdb/eggnog-mapper>

# Credits
>
Done by Alejandro Jimenez-Sanchez for his MSc dissertation.
>
The bioinformatics tools included in this repository were chosen based on the execution time and the quality of the results generated by each of them according to the bacterial genomes we worked with. Different programs were evaluated, mainly assemblers and annotators, which are not listed. It is recommended to test alternatives to those proposed here as each project is unique. 
