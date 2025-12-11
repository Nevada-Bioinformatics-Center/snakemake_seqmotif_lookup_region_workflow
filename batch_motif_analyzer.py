import argparse
import sys
import os
import glob
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

def load_motifs_from_tsv(tsv_path):
    """
    Loads sequence motifs from a two-column TSV file. 
    The first column is assumed to be the motif name, the second is the sequence.
    """
    if not os.path.exists(tsv_path):
        print(f"Error: Motif TSV file not found at {tsv_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading motifs from: {tsv_path}...")
    
    try:
        # Read the TSV file using tab delimiter
        # The header is assumed to be present and used for column names
        motif_df = pd.read_csv(tsv_path, sep='\t')
        
        # Ensure the DataFrame has at least two columns
        if motif_df.shape[1] < 2:
            print("Error: Motif TSV must contain at least two columns (Name and Sequence).", file=sys.stderr)
            sys.exit(1)

        # Use the first two columns for Name and Sequence
        motif_df = motif_df.iloc[:, 0:2].copy()
        motif_df.columns = ['Name', 'Sequence']
        
        # Convert DataFrame rows to a dictionary {Name: Sequence}
        #motifs = pd.Series(motif_df.Sequence.values, index=motif_df.Name).to_dict()
        motifs = pd.Series(motif_df.Sequence.astype(str).str.strip().values, index=motif_df.Name).to_dict()

        
        if not motifs:
            print("Error: No valid motifs were loaded from the TSV file.", file=sys.stderr)
            sys.exit(1)
            
        print(f"Successfully loaded {len(motifs)} motifs.")
        return motifs
        
    except Exception as e:
        print(f"Error reading or parsing motif TSV file: {e}", file=sys.stderr)
        sys.exit(1)


def analyze_fasta_for_motifs(fasta_path, motifs):
    """
    Analyzes a single FASTA file to find and count occurrences of defined motifs.
    Checks BOTH the forward sequence and its reverse complement.
    """
    total_reads = 0
    match_counts = defaultdict(int)
    
    # Pre-calculate Reverse Complements for efficiency
    # We store tuples of (Forward, ReverseComplement) for each motif name
    motif_variants = {}
    for name, seq in motifs.items():
        # Ensure input is string and uppercase
        fwd = str(seq).upper()
        # Calculate reverse complement using Biopython
        rc = str(Seq(fwd).reverse_complement())
        motif_variants[name] = (fwd, rc)
    
    try:
        # 1. Parse the FASTA file
        for record in SeqIO.parse(fasta_path, "fasta"):
            total_reads += 1
            # Convert sequence to string and uppercase for consistent searching
            read_seq = str(record.seq).upper()  

            # 2. Check for each motif in the read sequence
            for name, (fwd, rc) in motif_variants.items():
                # Check BOTH orientations
                if fwd in read_seq or rc in read_seq:
                    match_counts[name] += 1
            
    except Exception as e:
        print(f"Error processing {fasta_path}: {e}", file=sys.stderr)
        return None

    # 3. Compile results for this single sample
    sample_id = os.path.splitext(os.path.basename(fasta_path))[0]
    results = {'Sample_ID': sample_id, 'Total_Reads': total_reads}

    if total_reads == 0:
        print(f"Warning: {sample_id} contained 0 reads. Skipping motif calculations.")
        for name in motifs.keys():
            results[f'Count_{name}'] = 0
    else:
        for name in motifs.keys():
            count = match_counts[name]
            # Add counts to the result dictionary
            results[f'Count_{name}'] = count

    return results

def main():
    parser = argparse.ArgumentParser(description="Batch analyzer for genetic motifs.")
    
    parser.add_argument('-i', '--input', type=str, required=True, help='Directory containing FASTA files')
    parser.add_argument('-m', '--motifs', type=str, required=True, help='TSV file containing motif sequences')
    parser.add_argument('-o', '--output', type=str, default='motif_analysis_matrix.tsv', help='Output TSV matrix')

    args = parser.parse_args()
    
    motifs = load_motifs_from_tsv(args.motifs)
    motif_names = list(motifs.keys()) 
    
    search_path_fasta = os.path.join(args.input, '*.fasta')
    search_path_fa = os.path.join(args.input, '*.fa')
    fasta_files = glob.glob(search_path_fasta)
    fasta_files.extend(glob.glob(search_path_fa))

    if not fasta_files:
        print(f"Error: No .fasta or .fa files found in {args.input}", file=sys.stderr)
        sys.exit(1)
        
    print(f"Found {len(fasta_files)} FASTA files for analysis...")
    
    all_sample_data = []

    for i, fasta_file in enumerate(fasta_files):
        sample_results = analyze_fasta_for_motifs(fasta_file, motifs)
        if sample_results:
            all_sample_data.append(sample_results)

    if not all_sample_data:
        print("No samples were successfully analyzed.", file=sys.stderr)
        sys.exit(1)

    df = pd.DataFrame(all_sample_data)
    df = df.sort_values(by='Sample_ID', ascending=True).reset_index(drop=True)

    # Prepare columns for export
    cols = ['Sample_ID', 'Total_Reads']
    for name in motifs.keys():
        cols.append(f'Count_{name}')
    
    # Save Matrix
    df_final = df[cols]
    df_final.to_csv(args.output, index=False, sep='\t')
    print(f"Results saved to matrix: {os.path.abspath(args.output)}")


if __name__ == "__main__":
    try:
        import pandas as pd
        from Bio import SeqIO
        from Bio.Seq import Seq 
    except ImportError:
        print("Missing libraries. Install: pip install pandas biopython", file=sys.stderr)
        sys.exit(1)

    main()
