# snakemake_seqmotif_lookup_region_workflow

A snakemake workflow for extracting reads from a BAM file aligned to a certain region. Using a motifs.txt file, then per sample count the number of reads with certain motifs.

## Requirements

To run this pipeline, you will need [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) installed. Additionally, you will either need to install the required programs and packages on your system, or will need to run Snakemake using a container management system like Singularity, Apptainer or Docker.

### Programs Used

* samtools
* python
  * biopython
  * pandas
  * numpy
  * plotly


## Usage

#### Step 1: Install workflow

clone this workflow to your local computer

#### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`

#### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally (if you have all the software installed) via

    snakemake --cores $N

Execute the workflow locally (if you do not have all the software installed) via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-singularity --cluster qsub --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.
