# RNAseq sin genoma de referencia

Este pipeline realiza un análisis transcriptómico completo a partir de lecturas crudas en formato FASTQ. Está diseñado para ser reproducible, modular y ejecutarse tanto en local como en entornos HPC con Slurm.

## 1. Requisitos e instalación

Este proyecto usa [Snakemake](https://snakemake.readthedocs.io) y entornos `conda` para manejar dependencias. Si no tienes Snakemake instalado:

### Instalar Miniconda (si no lo tienes)
```bash
# En Linux
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```bash
### Crear entorno con snakemake

```bash
conda create -n snakemake_env -c conda-forge -c bioconda snakemake
conda activate snakemake_env
```

### Descarga de proyecto

Primero, crea un directorio para tu proyecto (si aún no lo has hecho) y clona el repositorio desde GitHub dentro de él:

```bash
mkdir -p mi_proyecto
cd mi_proyecto
git clone https://github.com/trippv/2025_snakemake_denovo.git
cd 2025_snakemake_denovo
```


## 2. Proceso
El flujo de trabajo automatiza las siguientes etapas:

1. Control de calidad y recorte: fastp

2. Construcción de ensamble de novo con Trinity

3. Estadisticos del ensamble con TrinityStats.pl

4. Control de calidad del ensamble con BUSCO

5. Control de calidad del ensamble con BOWTIE2

6. Cuantificación con Kallisto y generación de matrices


7. Generación de metadata y análisis PCA: script en Rmarkdown

8. Todos los pasos se ejecutan dentro de entornos conda definidos en la carpeta envs/.




## Estructura de archivos
```bash
project/
├── config/
│   ├── config.yaml        # Rutas de genoma y archivo de muestras
│   └── samples.tsv        # Información de muestras
├── data/                  # Carpeta de trabajo
├── testdata/              # (Ejemplo) Datos de prueba
├── results                # Salida del análisis
│   ├── busco
│   ├── Exp_Diff
│   ├── fastp
│   ├── kallisto
│   ├── mapping
│   ├── quant
│   ├── rnaquast
│   ├── summary_qc
│   ├── transdecoder
│   └── trinity               
├── rules/                 # Reglas individuales de Snakemake
├── envs/                  # Entornos conda por herramienta
├── scripts/               # Scripts auxiliares (e.g., Deseq2)
├── logs                   # Registro de sterr
└── Snakefile              # Punto de entrada del workflow
```

# Definir objetos

## lecturas

Crea el archivo samples.tsv dentro de `config/`. Ejemplo de contenido 


### Ejemplo de `samples.tsv`

| sample_id | group | fastq1                               | fastq2                           |extension    | include |
|-----------|--------|--------------------------------------|---------------------------------|-----|---------|
| batch1    | batch | testdata/raw/batch1_chrI_1.fastq     | testdata/raw/batch1_chrI_2.fastq |fastq    | 1       |
| chem1     | chem  | testdata/raw/chem1_chrI_1.fastq      | testdata/raw/chem1_chrI_2.fastq   |fastq   | 1       |


>include: define si la muestra se procesa (1) o se omite (0)

> extension: define la extension del archivo (.fastq, .fq, .fastq.gz, .fq.gz)

> Cada columna tiene que estar separada por tabulador (\t)

## Archivo de configuración

Archivo de configuración ubicado en `config/config.yaml`. Ejemplo

```bash
samples_file: "config/samples.tsv"

```

## Ejecutar el pipeline localmente

```bash
snakemake --use-conda --cores 8

```

## Ejecutar desde un administrador SLURM

Cuando se tienen limitantes en el acceso a multiples tareas en un administrador, la forma mas fácil de ejecutar snakemake es solicitando un solo nodo con el máximo de capacidad permitida. 

Un script se encuentra disponible en `scripts/run_snakemake_slurm_alt.sh`

```
#!/bin/bash
#SBATCH --job-name=RnaseqSnake
#SBATCH --output=rnaseq_snake_%j.log
#SBATCH --error=rnaseq_snake_%j.err
#SBATCH -p cicese
#SBATCH --ntasks-per-node=24
#SBATCH --mem=100GB
#SBATCH -t 6-00:00:00  # Cambia esto según lo que necesites

# Cargar Conda/Mamba
source /LUSTRE/apps/Miniforge/2024/miniforge3/etc/profile.d/conda.sh


conda activate snakemake  # cambia al entorno que uses

cd $SLURM_SUBMIT_DIR


# Ejecutar Snakemake en modo local con 24 núcleos
snakemake --cores 24 --use-conda


```