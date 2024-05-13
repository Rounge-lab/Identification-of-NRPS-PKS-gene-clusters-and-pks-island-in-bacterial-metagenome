# Pipeline for BGC Detection and Gene Annotation (not executed)

This pipeline is developed to detect biosynthetic gene clusters (BGCs) within metagenome assembled genomes (MAGs), and subsequently annotate the genes within the detected gene clusters. The BGC detection process is performed by the tools [antiSMASH](https://github.com/antismash/antismash) and [SanntiS](https://github.com/Finn-Lab/SanntiS), and the gene annotation process is performed by the tool [DRAM](https://github.com/WrightonLabCSU/DRAM).


## Installation

The necessary software and dependencies can be installed using the provided YAML files. The YAML files for the Conda environments are available [here](https://github.com/Rounge-lab/Identification-of-pks-positive-bacterial-genomes-in-CRCbiome/tree/main/pipeline1/envs).
```
conda env create -f [environment].yaml
```
- Replace "[environment]" with the actual file name of the YAML file.

The pipeline relies on certain Python packages and InterProscan for its execution, and these dependencies are loaded using module systems.


## Usage

```
snakemake -j [number_of_cores] --use-conda --conda-prefix [dir_conda_envs]
```
- Replace "[number_of_cores]" with the desired number of cores to be used by Snakemake for parallel execution.
- Replace "[directory_conda_environments]" with the actual path of your Conda environments directory.


## Configuration

The pipeline was designed to operate within [TSD](https://www.uio.no/english/services/it/research/sensitive-data/index.html), a secure platform for storage and analysis of sensitive data at the University of Oslo. For security reasons, sensitive data and their respective storage paths are not disclosed here. 


## Workflow Overview

The pipeline comprises a total of 10 rules, including 1 "rule all", 3 rules for the antiSMASH detection process and the subsequent DRAM annotation, and the remaining 6 rules for the SanntiS detection process and the subsequent DRAM annotation. Two custom Python scripts are incorporated for file format conversion and file modification. The Python scripts are available [here](https://github.com/Rounge-lab/Identification-of-pks-positive-bacterial-genomes-in-CRCbiome/tree/main/pipeline1/scripts).

![pipeline](https://github.com/Rounge-lab/Identification-of-pks-positive-bacterial-genomes-in-CRCbiome/blob/main/figures/pipeline1.png)
