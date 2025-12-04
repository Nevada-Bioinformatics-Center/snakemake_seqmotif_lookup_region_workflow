configfile: "config.yaml"

import os
import glob
import sys
import collections

SAMTOOLS_IMAGE = "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1" 
PYTHON_IMAGE = "docker://community.wave.seqera.io/library/pip_biopython_pandas_plotly:6c94c410be71f235"

SAMPLE_TO_BAM_FILE = collections.OrderedDict()
# Use glob to find all BAM files in the input directory
BAM_FILES = glob.glob(os.path.join(config["INPUT_DIR"], "*.bam")) 

SUFFIXES = [".sortedByCoord.out.bam", ".sorted.bam"]
DELIMITER_NEW = config['DELIMITER']
REPLICATE_SUFFIX = ".rep1" 

for f in BAM_FILES:
    basename = os.path.basename(f)
    
    clean_name = basename
    for suffix in SUFFIXES:
        if clean_name.endswith(suffix):
            clean_name = clean_name.replace(suffix, "")
            break
    
    sample_id = None

    if DELIMITER_NEW in clean_name:
        sample_id = clean_name.split(DELIMITER_NEW)[0]
    else:
        sample_id = clean_name
        if sample_id.endswith(REPLICATE_SUFFIX):
            sample_id = sample_id.replace(REPLICATE_SUFFIX, "")
        
    if sample_id:
        if sample_id in SAMPLE_TO_BAM_FILE:
            print(f"Error: Duplicate sample ID '{sample_id}' found for files {SAMPLE_TO_BAM_FILE[sample_id]} and {basename}. Exiting.", file=sys.stderr)
            sys.exit(1)
            
        SAMPLE_TO_BAM_FILE[sample_id] = basename
    else:
        print(f"Warning: Could not determine sample ID for file {basename}. Skipping.", file=sys.stderr)


SAMPLES = list(SAMPLE_TO_BAM_FILE.keys())

if not SAMPLES:
    print(f"Error: No BAM files found in {config['INPUT_DIR']}. Exiting.", file=sys.stderr)
    sys.exit(1)


def get_full_bam_path(wildcards):
    """Returns the full path for a BAM file given the short {sample} ID."""
    full_filename = SAMPLE_TO_BAM_FILE[wildcards.sample]
    return os.path.join(config["INPUT_DIR"], full_filename)

def get_full_bam_index_path(wildcards):
    """
    Calculates the required index path: [FULL_BAM_PATH].bai.
    """
    return get_full_bam_path(wildcards) + ".bai"


# --- Workflow Rules ---

rule all:
    input:
        os.path.join(config["WORK_DIR"], config["OUTPUT_MATRIX"])

rule index_bam:
    input:
        os.path.join(config["INPUT_DIR"], "{filename}.bam")
    output:
        os.path.join(config["INPUT_DIR"], "{filename}.bam.bai")
    #conda: "envs/samtools.yaml"
    container: SAMTOOLS_IMAGE 
    threads: 2
    resources: time_min=220, mem_mb=8000, cpus=2
    shell:
        """
        samtools index {input} -o {output}
        """

rule extract_and_convert:
    input:
        bam = get_full_bam_path,
        index = get_full_bam_index_path 
    output:
        fasta = os.path.join(config["WORK_DIR"], "fasta_reads", "{sample}.fasta")
    params:
        region = config["TARGET_REGION"]
    log:
        os.path.join(config["WORK_DIR"], "logs", "{sample}.extract.log")
    #conda: "envs/samtools.yaml"
    container: SAMTOOLS_IMAGE 
    threads: 2
    resources: time_min=220, mem_mb=8000, cpus=2
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        
        # samtools view finds the index automatically because it is right next to {input.bam}
        samtools view -h {input.bam} {params.region} | 
        samtools fasta - > {output.fasta} 2> {log}
        """

rule analyze_all:
    input:
        fasta_files = [
            os.path.join(config["WORK_DIR"], "fasta_reads", f"{sample}.fasta") 
            for sample in SAMPLES
        ],
        motifs = config["MOTIF_TSV"]
    output:
        matrix = os.path.join(config["WORK_DIR"], config["OUTPUT_MATRIX"]),
    params:
        input_dir = os.path.join(config["WORK_DIR"], "fasta_reads")
    log:
        os.path.join(config["WORK_DIR"], "logs", "final_analysis.log")
    #conda: "envs/python_analysis.yaml"
    container: PYTHON_IMAGE
    threads: 2
    resources: time_min=220, mem_mb=8000, cpus=2
    shell:
        """
        mkdir -p $(dirname {output.matrix})
        
        # No need to install anything; the image has everything ready.
        python batch_motif_analyzer.py \\
            -i {params.input_dir} \\
            -m {input.motifs} \\
            -o {output.matrix} 2> {log}
        """
