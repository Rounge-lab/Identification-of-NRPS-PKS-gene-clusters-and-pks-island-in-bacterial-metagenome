# MAG-based identification of BGCs

This pipeline is developed to detect biosynthetic gene clusters (BGCs) within metagenome assembled genomes (MAGs), annotate the genes within the detected gene clusters, and subsequently extract the resulting data for further analysis. The BGC detection process is performed by the tool [antiSMASH](https://github.com/antismash/antismash), and the gene annotation process is performed by the tool [DRAM](https://github.com/WrightonLabCSU/DRAM).


## Installation

The necessary software and dependencies can be installed using the provided YAML files. The YAML files for the Conda environments are available [here](https://github.com/Rounge-lab/Identification-of-NRPS-PKS-gene-clusters-and-pks-island-in-bacterial-metagenome/tree/main/pipeline2/envs).
```
conda env create -f [environment].yaml
```
- Replace "[environment]" with the actual file name of the YAML file.

The pipeline relies on certain Python packages for its execution, and these dependencies are loaded using module systems.


## Usage

The pipeline offers two execution options: locally or via SLURM.

**Local execution**
```
snakemake -j [number_of_cores] --use-conda --conda-prefix [dir_conda_envs]
```
- Replace "[number_of_cores]" with the desired number of cores to be used by Snakemake for parallel execution.
- Replace "[directory_conda_environments]" with the actual path of your Conda environments directory.

**SLURM execution**

```<div style="overflow-x: auto; white-space: nowrap;">
snakemake -j [number_of_cores] --use-conda --conda-prefix [dir_conda_envs] --profile config/cluster/ --latency-wait [number of sec]
```
- Replace "[number_of_cores]" with the desired number of cores to be used by Snakemake for parallel execution.
- Replace "[directory_conda_environments]" with the actual path of your Conda environments directory.
- The folder [config/cluster/](https://github.com/Rounge-lab/Identification-of-NRPS-PKS-gene-clusters-and-pks-island-in-bacterial-metagenome/tree/main/pipeline2/config/cluster) contains the profile configuration files that customize Snakemake's behavior for specific execution environments.
- Replace "[number of sec]" with the desired latency wait time in seconds, which is the time Snakemake waits for file system changes before considering a job ready to be executed.


## Configuration

The pipeline was designed to operate within [TSD](https://www.uio.no/english/services/it/research/sensitive-data/index.html), a secure platform for storage and analysis of sensitive data at the University of Oslo. For security reasons, sensitive data and their respective storage paths are not disclosed here. 


## Workflow Overview

The pipeline comprises a total of 12 rules, including 1 "rule all", 5 rules for the core detection and annotation process, and the remaining 6 rules for extracting and concatenating output data for further data analysis. Eight custom Python scripts are incorporated for file format conversion, data extraction, and file concatenation. The Python scripts are available [here](https://github.com/Rounge-lab/Identification-of-NRPS-PKS-gene-clusters-and-pks-island-in-bacterial-metagenome/tree/main/pipeline2/scripts).

![pipeline](https://github.com/Rounge-lab/Identification-of-NRPS-PKS-gene-clusters-and-pks-island-in-bacterial-metagenome/blob/main/flowchart/pipeline2.png)
