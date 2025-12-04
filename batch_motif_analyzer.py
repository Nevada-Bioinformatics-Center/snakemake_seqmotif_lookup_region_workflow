import argparse
import sys
import os
import glob
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
#import plotly.express as px
#import plotly.graph_objects as go


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
        motifs = pd.Series(motif_df.Sequence.values, index=motif_df.Name).to_dict()
        
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
    """
    total_reads = 0
    match_counts = defaultdict(int)
    
    try:
        # Parse the FASTA file
        for record in SeqIO.parse(fasta_path, "fasta"):
            total_reads += 1
            # Convert sequence to string and uppercase for consistent searching
            read_seq = str(record.seq).upper()  

            #  Check for each motif in the read sequence
            for name, seq in motifs.items():
                # Perform a direct, exact-match search
                if seq in read_seq:
                    match_counts[name] += 1
            
    except Exception as e:
        print(f"Error processing {fasta_path}: {e}", file=sys.stderr)
        return None

    # Compile results for this single sample
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

#def assign_category(row, motif_names):
#    """
#    Assigns a category based on the counts of the motifs.
#    Rules:
#    1. If all count columns >= 1 -> "All_motif_hits"
#    2. If first category >= 1 AND all others == 0 -> "WTonly_motif_hits" (Assuming first is WT)
#    3. If first category == 0 AND at least one other >= 1 -> "NoWT_motif_hits"
#    4. Else -> "Other_motif_hits"
#    """
#    # Extract just the count values for the motifs
#    # Use .get() to safely handle potential missing columns, defaulting to 0
#    counts = [row.get(f'Count_{name}', 0) for name in motif_names]
#    
#    if not counts:
#        return "Other_motif_hits"
#
#    first_count = counts[0]
#    other_counts = counts[1:]
#
#    # Rule 1: All columns have at least a count of 1
#    if all(c >= 1 for c in counts):
#        return "All_motif_hits"
#    
#    # Rule 2: All counts but the first category are 0 (and first is present)
#    if first_count >= 1 and all(c == 0 for c in other_counts):
#        return "WTonly_motif_hits"
#    
#    # Rule 3: WT is 0, but we have hits in other columns
#    if first_count == 0 and any(c >= 1 for c in other_counts):
#        return "NoWT_motif_hits"
#
#    # Rule 4: Everything else (e.g. Mixed hits but not ALL, or No hits at all)
#    return "Other_motif_hits"

#def generate_plot(df, motif_names, output_path):
#    """Generates an interactive Plotly scatter plot using Counts."""
#    
#    # Calculate axes for plotting
#    # X-Axis: Count of the First Motif (WT)
#    x_axis_col = f'Count_{motif_names[0]}'
#    
#    # Y-Axis: Sum of Counts of ALL other motifs
#    y_vals = df[[f'Count_{name}' for name in motif_names[1:]]].sum(axis=1)
#    df['Count_Other_Motifs'] = y_vals
#
#    print("Generating Plotly figure...")
#    
#    # Define the list of columns we want to see on hover
#    hover_cols = ['Sample_ID', 'Total_Reads']
#    # Add individual motif counts
#    for name in motif_names:
#        hover_cols.append(f'Count_{name}')
#    
#    MAX_BUBBLE_SIZE = 15
#
#    # Create standard Plotly Express scatter plot
#    fig = px.scatter(
#        df,
#        x=x_axis_col,
#        y='Count_Other_Motifs',
#        color='WT_check',
#        hover_data=hover_cols,
#        title='Genetic Motif Distribution by Sample (Counts)',
#        labels={
#            x_axis_col: "WT Motif Hits",
#            "Count_Other_Motifs": "All Mutant Motif Hits",
#            "WT_check": "Motif Type",
#            "Total_Reads": "Total Reads"
#        },
#        color_discrete_map={
#            "All_motif_hits": "purple", 
#            "WTonly_motif_hits": "blue",
#            "NoWT_motif_hits": "red",
#            "Other_motif_hits": "orange"
#        },
#        size='Total_Reads', 
#        size_max=MAX_BUBBLE_SIZE
#    )
#    
#    max_reads = df['Total_Reads'].max()
#    min_reads = df['Total_Reads'].min()
#    
#    if pd.isna(max_reads) or max_reads == 0:
#        max_reads = 100
#        min_reads = 10
#
#    fig.add_trace(go.Scatter(
#        x=[None], y=[None],
#        mode='markers',
#        marker=dict(size=2, color='grey'), 
#        name=f'{int(min_reads)} Reads (Min)',
#        showlegend=True,
#        legendgroup='size'
#    ))
#
#    fig.add_trace(go.Scatter(
#        x=[None], y=[None],
#        mode='markers',
#        marker=dict(size=MAX_BUBBLE_SIZE, color='grey'), 
#        name=f'{int(max_reads)} Reads (Max)',
#        showlegend=True,
#        legendgroup='size'
#    ))
#
#    fig.update_layout(
#        legend=dict(itemsizing='trace', title_text='Legend')
#    )
#
#    fig.write_html(output_path)
#    print(f"Plot saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Batch analyzer for genetic motifs.")
    
    parser.add_argument('-i', '--input', type=str, required=True, help='Directory containing FASTA files')
    parser.add_argument('-m', '--motifs', type=str, required=True, help='TSV file containing motif sequences')
    parser.add_argument('-o', '--output', type=str, default='motif_analysis_matrix.tsv', help='Output TSV matrix')
    #parser.add_argument('--plot', type=str, default=None, help='Path to save interactive HTML plot')

    args = parser.parse_args()
    
    motifs = load_motifs_from_tsv(args.motifs)
    motif_names = list(motifs.keys()) # Keep order for categorization
    
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

    #df['WT_check'] = df.apply(lambda row: assign_category(row, motif_names), axis=1)

    cols = ['Sample_ID', 'Total_Reads']
    for name in motifs.keys():
        cols.append(f'Count_{name}')
    
    # Add the new category column to the export list
    #cols.append('WT_check')
    
    # Save Matrix
    df_final = df[cols]
    df_final.to_csv(args.output, index=False, sep='\t')
    print(f"Results saved to matrix: {os.path.abspath(args.output)}")

    # Generate Plot if requested 
    #if args.plot:
    #    generate_plot(df_final, motif_names, args.plot)


if __name__ == "__main__":
    try:
        import pandas as pd
        from Bio import SeqIO
        import plotly
    except ImportError:
        print("Missing libraries. Install: pip install pandas biopython plotly", file=sys.stderr)
        sys.exit(1)

    main()
