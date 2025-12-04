# snakemake_seqmotif_lookup_region_workflow

A snakemake workflow for extracting reads from a BAM file aligned to a certain region. Using a motifs.txt file, calculate the number of reads with defined motifs per-sample.

## Requirements

To run this pipeline, you will need [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) installed. Additionally, you will either need to install the required programs and packages on your system or use the `--use-conda` flag to install the conda environemtns or use the `--use-singularity` flag to run Snakemake using a container management system like Singularity or Apptainer (typically only on HPCs).

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

## batch_motif_analyzer.py

This, [`batch_motif_analyzer.py`](batch_motif_analyzer.py), is the main analysis script. The snakemake rules beforehand takes in a directory of alignment BAM files and creates .fasta files for the reads within extracted region.

The batch_motif_analyzer.py takes in as input a directory location of extracted reads (`--input`) and a motifs file in TSV format (`--motifs`) to run the analysis. The motifs file must have a header line where the first column is the name of the motif and the second column is the sequence of the motif to search. Additionally, the analysis assumes the first entry is the WT motif which will be labeled as such in the output `Plotly` HTML output report. See example file: [`motifs.txt`](motifs.txt)

For output, use the `--output` option to name the output matrix of counts for each Motif and sample. Optionally, define a `--plot` HTML file to save an interactive output of the WT readcounts x Other readcounts.

#### Example Usage

`python batch_motif_analyzer.py --input reads_fasta/ --motifs motifs.txt --output my_counts_matrix.tsv --plot motif_analysis.html`
