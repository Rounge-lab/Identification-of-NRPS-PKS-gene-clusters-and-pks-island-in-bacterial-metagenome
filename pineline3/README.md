# Read-based identification of pks island and E. coli

This pipeline is developed to align sequencing libraries agains reference genomes of interest, generate mapping statitics, and extract the resulting data for further analysis. These processes utilized the tools [Bowtie2](https://github.com/BenLangmead/bowtie2) and [SAMtools](https://github.com/samtools/samtools). 


## Installation

The necessary software and dependencies can be installed using the provided YAML files. The YAML files for the Conda environments are available [here](https://github.com/Rounge-lab/Identification-of-pks-positive-bacterial-genomes-in-CRCbiome/tree/main/pineline3/envs).
```
conda env create -f [environment].yaml
```
- Replace "[environment]" with the actual file name of the YAML file.

The pipeline relies on certain Python packages for its execution, and these dependencies are loaded using module systems.


## Usage

```
snakemake -j [number_of_cores] --use-conda --conda-prefix [dir_conda_envs]
```
- Replace "[number_of_cores]" with the desired number of cores to be used by Snakemake for parallel execution.
- Replace "[directory_conda_environments]" with the actual path of your Conda environments directory.


## Configuration

The pipeline was designed to operate within [TSD](https://www.uio.no/english/services/it/research/sensitive-data/index.html), a secure platform for storage and analysis of sensitive data at the University of Oslo. For security reasons, sensitive data and their respective storage paths are not disclosed here. 


## Workflow Overview

The pipeline comprises a total of 11 rules, including 1 "rule all", 2 rules for the core mapping process utilizing Bowtie2 and SAMtools, 2 for obtaining mapping statistics utilizing SAMtools, and the remaining 6 rules for organizing the statistical data for further analysis. Three custom Python scripts are incorporated for data extraction and file concatenation, and two BED files are required for generating mapping statistics. The Python scripts and BED files are available [here](https://github.com/Rounge-lab/Identification-of-pks-positive-bacterial-genomes-in-CRCbiome/tree/main/pineline3/scripts).

![pipeline](https://github.com/Rounge-lab/Identification-of-pks-positive-bacterial-genomes-in-CRCbiome/blob/main/figures/pipeline3.png)
