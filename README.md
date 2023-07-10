> **TAIFA is a workflow to analyse whole genome sequencing data from bacteria**
>
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

# Introduction

# Workflow summary
- Sequencing quality control and trimming
- De novo genome assembly and evaluation 
- Alignment 
- Annotation
- Pathogenicity prediction 
  
# Workflow
## 1. Sequencing quality control and trimming
1.1. First FastQC
```{bash}
/path/to/fastqc path/to/*.fastq.gz -o "$carpeta"/Alignment_1/Fastq/fastqc
```
1.2. Trimmomatic
```{bash}
# Ruta de la carpeta principal
java -jar /data/home/alejandro_jimenez/trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 MG_S1_L001_R1_001.fastq.gz MG_S1_L001_R2_001.fastq.gz /home/alejandro_jimenez/seq/mg/Alignment_1/Fastq/fastqc/trim_MG_S1_L001_R1_001.fastq.gz Undetermined_S0_L001_R1_001.fastq.gz /home/alejandro_jimenez/seq/mg/Alignment_1/Fastq/fastqc/trim_MG_S1_L001_R2_001.fastq.gz Undetermined_S0_L001_R2_001.fastq.gz CROP:150 HEADCROP:1 AVGQUAL:20 SLIDINGWINDOW:4:20
```
1.3. Second FastQC
```{bash}
/path/to/fastqc path/to/*.fastq.gz -o "$carpeta"/Alignment_1/Fastq/fastqc
```
## 2. De novo genome assembly and evaluation 
2.1. A5-miseq
>
Description:
>
Link: https://sourceforge.net/p/ngopt/wiki/A5PipelineREADME/  
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
          nohup time -p /data/home/alejandro_jimenez/a5/bin/a5_pipeline.pl --threads=4 "$archivo_r1" "$archivo_r2" a5_output &
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

2.2. QUAST
>
Description:
>
Link: https://quast.sourceforge.net/docs/manual.html
```{bash}
conda activate quast_env
cd /home/alejandro_jimenez/seq/mg/Alignment_1/Fastq/fastqc
quast.py a5_output.contigs.fasta -o quast_a5_output -l a5
````

2.3. BUSCO
>
Description: 
>
Link: https://busco.ezlab.org/
```{bash}
module load Anaconda3/4.4.0
module load augustus
run_BUSCO.py -c 4 -i a5_output.contigs.fasta -l /home/alejandro_jimenez/bacteria_odb9 -o busco_a5 -m geno 
```

## 3. Mapping
3.1. BWA, SAMtools
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

3.2. IGV
>
Description:
>
Link:

3.3. Qualimap

## 4. Annotation
4.1. Prokka
>
Description:
>
Link: https://github.com/tseemann/prokka/blob/master/README.md
```{bash}
conda activate prokka_env
cd /home/alejandro_jimenez/seq/mg/Alignment_1/Fastq/fastqc/
prokka a5_output.contigs.fasta --cpus 6 --outdir prokka_output
```

4.2. Functional annotation
>
Description:
>
Link:
http://eggnog-mapper.embl.de/

## 5. Virulence factors prediction
>
Description:
>
Link:
```{bash
conda activate diamond_env
diamond blastp --db /home/alejandro_jimenez/vfdb -q SA.fasta -o vfdb_results_SA.tsv --query-cover 97 --id 50
```

# Credits
