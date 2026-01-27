#!/usr/bin/env python3
"""
MS-RT Benchmark for Hyperspec
Calculate N10 and Purity (MS-RT validation) metrics for GPU clustering results.
"""

import os
import sys
import pandas as pd
import numpy as np
import networkx as nx
from collections import Counter, defaultdict
from itertools import combinations
from multiprocessing import Pool
from functools import partial
from tqdm import tqdm
import gc
import argparse
import subprocess
from pathlib import Path
import re

# Method dictionary for falcon format
method_dic = {
    'falcon': {
        'filename': 'filename',
        'scan': 'scan',
        'mass': 'precursor_mz',
        'rt_time': 'retention_time'
    }
}


def extract_dataset_name(file_path):
    """Extract dataset name from file path"""
    # Try to extract MSV number or dataset identifier
    path_str = str(file_path)
    
    # Look for MSV number pattern
    msv_match = re.search(r'MSV\d+', path_str)
    if msv_match:
        return msv_match.group(0)
    
    # Look for other patterns
    # Extract from filename like "MSV000081981_clustering_results_1k.csv.parquet"
    filename = Path(file_path).stem  # Remove .parquet
    filename = Path(filename).stem    # Remove .csv if present
    
    # Extract dataset prefix
    match = re.match(r'([^_]+)', filename)
    if match:
        return match.group(1)
    
    # Fallback: use parent directory name or generic name
    parent = Path(file_path).parent.name
    if parent:
        return parent
    
    return "dataset"


def convert_parquet_to_tsv(input_file):
    """
    Convert parquet file to TSV format using convert_to_tsv.py
    
    Args:
        input_file: Path to parquet file
        
    Returns:
        Path to converted TSV file
    """
    input_path = Path(input_file)
    
    # Auto-generate output filename using input file prefix + .tsv
    # e.g., orbitrap_clustering_results_1k.csv.parquet -> orbitrap_clustering_results_1k.csv.tsv
    output_dir = input_path.parent
    # Remove .parquet or .csv.parquet extension and add .tsv
    stem = input_path.stem  # Gets filename without extension
    if stem.endswith('.csv'):
        # Handle .csv.parquet case
        stem = stem[:-4]  # Remove .csv
    output_path = output_dir / f"{stem}.tsv"
    
    # Check if converted file already exists
    if output_path.exists():
        print(f"Converted TSV file already exists: {output_path}")
        print("Using existing file. Delete it if you want to regenerate.")
        return str(output_path)
    
    print(f"Detected parquet file: {input_path}")
    print(f"Converting to TSV format...")
    
    # Get the directory of this script
    script_dir = Path(__file__).parent
    convert_script = script_dir / "convert_to_tsv.py"
    
    if not convert_script.exists():
        # Fallback: try direct conversion using pandas
        print("convert_to_tsv.py not found, using direct pandas conversion...")
        try:
            df = pd.read_parquet(input_path)
            df.to_csv(output_path, sep='\t', index=False)
            print(f"✓ Conversion completed: {output_path}")
            return str(output_path)
        except Exception as e:
            print(f"Error: Failed to convert parquet file: {e}", file=sys.stderr)
            print("Please install pyarrow: pip install pyarrow", file=sys.stderr)
            sys.exit(1)
    
    # Use convert_to_tsv.py script
    try:
        result = subprocess.run(
            [sys.executable, str(convert_script), str(input_path), str(output_path)],
            capture_output=True,
            text=True,
            check=True
        )
        print(result.stdout)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        return str(output_path)
    except subprocess.CalledProcessError as e:
        print(f"Error: Failed to convert parquet file: {e}", file=sys.stderr)
        if e.stdout:
            print(e.stdout, file=sys.stderr)
        if e.stderr:
            print(e.stderr, file=sys.stderr)
        sys.exit(1)


def preprocess_input_file(input_file):
    """
    Preprocess input file: convert parquet to TSV if needed
    
    Args:
        input_file: Path to input file (TSV or parquet)
        
    Returns:
        Path to TSV file (converted if needed, original if already TSV)
    """
    input_path = Path(input_file)
    
    # Check if file is parquet format
    if input_path.suffix == '.parquet' or input_path.name.endswith('.csv.parquet'):
        print(f"\n{'='*80}")
        print("PREPROCESSING: Converting parquet to TSV format")
        print(f"{'='*80}")
        tsv_file = convert_parquet_to_tsv(input_file)
        print(f"{'='*80}\n")
        return tsv_file
    
    # Already TSV, return as is
    return input_file


def load_gpu_results(gpu_file):
    """
    Load GPU clustering results
    
    Returns:
        tuple: (filtered_dataframe, total_scans_in_file)
            - filtered_dataframe: DataFrame with cluster == -1 filtered out (for purity calculation)
            - total_scans_in_file: Total number of scans including unclustered (for N10 calculation)
    
    Note: scan == -1 indicates unclustered spectra (noise), which are filtered out
    for purity calculation but should be included in total scans for N10 calculation.
    """
    print(f"Loading GPU results from {gpu_file}...")
    df = pd.read_csv(gpu_file, sep='\t')
    
    # Save total scans before filtering (includes unclustered spectra with cluster == -1)
    total_scans_in_file = len(df)
    
    # GPU identifier format: 000011026_RA3_01_6002 (this is the base filename)
    # scan column is read directly from the file (may be -1)
    # For purity calculation, we'll use identifier as filename
    df['base_filename'] = df['identifier'].astype(str)
    df['filename'] = df['identifier'].astype(str) + '.mzML'  # Add .mzML for consistency
    
    # Use scan column directly from file
    # scan column is read as-is from the TSV file
    
    # Create spectrum identifier: identifier_precursor_mz_retention_time
    # This should be unique for each spectrum
    df['spectrum_id'] = (
        df['identifier'].astype(str) + '_' + 
        df['precursor_mz'].astype(str) + '_' + 
        df['retention_time'].astype(str)
    )
    
    # Filter out noise (cluster == -1) for purity calculation
    # These unclustered spectra are excluded from clustering metrics but should be
    # included in total scans for N10 calculation
    df_filtered = df[df['cluster'] != -1].copy()
    
    print(f"Loaded {len(df_filtered)} GPU scan assignments (clustered)")
    print(f"Total scans in file (including unclustered): {total_scans_in_file:,}")
    print(f"GPU clusters: {df_filtered['cluster'].nunique()}")
    
    return df_filtered, total_scans_in_file


def optimized_create_matching_network(cluster, method, rt_window=30.0, precursor_mz_window=0.01):
    """
    Create matching network for a cluster (optimized version)
    
    Args:
        cluster: DataFrame with cluster data
        method: Method name (e.g., 'falcon')
        rt_window: Retention time window in seconds (default: 30.0)
        precursor_mz_window: Precursor m/z window in Da (default: 0.01)
    """
    G = nx.Graph()

    # Precompute the node names and add them to the graph
    node_attrs = {
        f"{row[method_dic[method]['filename']]}_{row[method_dic[method]['scan']]}": {
            "filename": row[method_dic[method]['filename']]
        }
        for _, row in cluster.iterrows()
    }
    G.add_nodes_from(node_attrs.items())

    # Precompute mass and rt_time for efficient access
    specs = cluster.apply(lambda row: (
        f"{row[method_dic[method]['filename']]}_{row[method_dic[method]['scan']]}",
        row[method_dic[method]['mass']],
        row[method_dic[method]['rt_time']]
    ), axis=1).tolist()

    # Create edges based on conditions
    # RT tolerance: configurable (default: 30 seconds)
    # Precursor m/z tolerance: configurable (default: 0.01 Da)
    edges = [
        (spec1[0], spec2[0]) for spec1, spec2 in combinations(specs, 2)
        if spec1[0].split('_')[0] == spec2[0].split('_')[0]  # Ensure filenames are the same
           and abs(spec1[1] - spec2[1]) <= precursor_mz_window 
           and abs(spec1[2] - spec2[2]) <= rt_window
    ]
    G.add_edges_from(edges)

    return G


def calculate_max_component_per_file(G):
    """Calculate maximum component size for each file"""
    # Find all connected components in the graph
    components = nx.connected_components(G)

    # Initialize a dictionary to hold the maximum component size for each file
    max_component_sizes = defaultdict(int)

    # Iterate through each component
    for component in components:
        # Create a temporary dictionary to count the number of nodes per file in this component
        file_counts = defaultdict(int)

        # Count nodes per file in the current component
        for node in component:
            filename = G.nodes[node]['filename']
            file_counts[filename] += 1

        # Update the max component size for each file encountered in this component
        for filename, count in file_counts.items():
            if count > max_component_sizes[filename]:
                max_component_sizes[filename] = count

    # Ensure that files represented by single nodes are accounted for
    for node in G.nodes:
        filename = G.nodes[node]['filename']
        if filename not in max_component_sizes:
            max_component_sizes[filename] = 1
        else:
            # Ensure there's at least a count of 1 for each file
            max_component_sizes[filename] = max(max_component_sizes[filename], 1)

    return max_component_sizes


def process_cluster(cluster_id, cluster_data, method='falcon', rt_window=30.0, precursor_mz_window=0.01):
    """
    Process a single cluster to calculate purity (MS-RT method)
    
    Args:
        cluster_id: Cluster ID
        cluster_data: DataFrame with cluster data
        method: Method name (e.g., 'falcon')
        rt_window: Retention time window in seconds
        precursor_mz_window: Precursor m/z window in Da
    """
    if len(cluster_data) == 1:
        return (1.0, 1)  # Singleton cluster has purity 1
    
    G = optimized_create_matching_network(cluster_data, method, rt_window, precursor_mz_window)
    max_component_sizes = calculate_max_component_per_file(G)
    
    # Calculate the count of each filename in the cluster
    file_counts = Counter(cluster_data[method_dic[method]['filename']])
    
    # Calculate the fraction of the largest component for each file
    frequencies = []
    values = []
    for filename, count in file_counts.items():
        largest_component_size = max_component_sizes.get(filename, 1)  # Default to 1 if not found
        fraction = largest_component_size / count
        frequencies.append(count)
        values.append(fraction)
    
    # Calculate weighted average purity
    weighted_sum = sum(value * frequency for value, frequency in zip(values, frequencies))
    total_frequency = sum(frequencies)
    weighted_average = weighted_sum / total_frequency if total_frequency > 0 else 0.0
    
    cluster_size = len(cluster_data)
    return (weighted_average, cluster_size)


def falcon_purity(cluster_results, batch_size=10000, rt_window=30.0, precursor_mz_window=0.01):
    """
    Calculate purity for each cluster using network-based approach.
    Uses batching to avoid loading all clusters into memory at once.
    
    Args:
        cluster_results: DataFrame with cluster assignments
        batch_size: Number of clusters to process in each batch (default: 10000)
        rt_window: Retention time window in seconds (default: 30.0)
        precursor_mz_window: Precursor m/z window in Da (default: 0.01)
    
    Returns:
        purity_list: List of purity values for each cluster
        size_list: List of cluster sizes
    """
    method = 'falcon'
    purity_list = []
    size_list = []
    
    # Group by cluster (returns iterator, not list)
    clusters = cluster_results.groupby('cluster')
    
    # Get total number of clusters
    total_clusters = clusters.ngroups
    print(f"Processing {total_clusters:,} clusters in batches of {batch_size:,}...")
    print(f"Using RT window: {rt_window} seconds, Precursor m/z window: {precursor_mz_window} Da")
    
    # Batch processing to avoid memory issues
    batch_items = []
    processed_count = 0
    
    # Progress bar
    pbar = tqdm(total=total_clusters, desc="Calculating purity")
    
    # Process clusters in batches
    for cluster_id, cluster_data in clusters:
        batch_items.append((cluster_id, cluster_data))
        
        # Process batch when it reaches batch_size or at the end
        if len(batch_items) >= batch_size:
            # Process current batch with parallel processing
            with Pool() as pool:
                process_func = partial(process_cluster, method=method, 
                                     rt_window=rt_window, precursor_mz_window=precursor_mz_window)
                batch_results = pool.starmap(process_func, batch_items)
            
            # Collect results
            for purity, size in batch_results:
                purity_list.append(purity)
                size_list.append(size)
            
            processed_count += len(batch_items)
            pbar.update(len(batch_items))
            
            # Clear batch and force garbage collection
            batch_items = []
            del batch_results
            gc.collect()
    
    # Process remaining clusters (if any)
    if batch_items:
        with Pool() as pool:
            process_func = partial(process_cluster, method=method,
                                 rt_window=rt_window, precursor_mz_window=precursor_mz_window)
            batch_results = pool.starmap(process_func, batch_items)
        
        for purity, size in batch_results:
            purity_list.append(purity)
            size_list.append(size)
        
        pbar.update(len(batch_items))
        del batch_results
        gc.collect()
    
    pbar.close()
    print(f"Processed {len(purity_list):,} clusters successfully.")
    
    return purity_list, size_list


def calculate_n10(cluster_sizes, total_scans):
    """Calculate N10: cluster size at which 10% of total scans are covered"""
    # Sort cluster sizes in descending order
    sorted_sizes = sorted(cluster_sizes.values, reverse=True)
    
    # Calculate cumulative coverage
    cumulative = 0
    target_coverage = total_scans * 0.10  # 10% of total scans
    
    for size in sorted_sizes:
        cumulative += size
        if cumulative >= target_coverage:
            return size
    
    # If we never reach 10%, return the largest cluster size
    return sorted_sizes[0] if sorted_sizes else 0


def calculate_weighted_average_purity(purity_array, size_array):
    """Calculate weighted average purity"""
    if len(purity_array) == 0 or len(size_array) == 0:
        return 0.0
    
    total_weight = np.sum(size_array)
    if total_weight == 0:
        return 0.0
    
    weighted_sum = np.sum(purity_array * size_array)
    return weighted_sum / total_weight


def calculate_metrics(df, total_scans, rt_window=30.0, precursor_mz_window=0.01, batch_size=10000):
    """
    Calculate N10 and Purity metrics for a clustering result
    
    Args:
        df: DataFrame with cluster assignments
        total_scans: Total number of scans for N10 calculation
        rt_window: Retention time window in seconds (default: 30.0)
        precursor_mz_window: Precursor m/z window in Da (default: 0.01)
        batch_size: Batch size for processing (default: 10000)
    """
    print(f"\n{'='*80}")
    print(f"Calculating metrics for GPU clustering results")
    print(f"{'='*80}")
    
    # Calculate cluster sizes
    cluster_sizes = df.groupby('cluster').size()
    
    # Calculate N10
    print("Calculating N10...")
    n10 = calculate_n10(cluster_sizes, total_scans)
    print(f"N10 value: {n10:,}")
    
    # Calculate purity
    print("Calculating purity (this may take a while)...")
    purity_list, size_list = falcon_purity(df, batch_size=batch_size, 
                                          rt_window=rt_window, 
                                          precursor_mz_window=precursor_mz_window)
    
    # Convert to numpy arrays for easier calculation
    purity_array = np.array(purity_list)
    size_array = np.array(size_list)
    
    # Calculate weighted average purity
    weighted_avg_purity = calculate_weighted_average_purity(purity_array, size_array)
    print(f"Weighted Average Purity: {weighted_avg_purity:.6f}")
    
    # Calculate simple average purity
    simple_avg_purity = np.mean(purity_array)
    print(f"Simple Average Purity: {simple_avg_purity:.6f}")
    
    return {
        'n10': n10,
        'weighted_avg_purity': weighted_avg_purity,
        'simple_avg_purity': simple_avg_purity,
        'total_clusters': len(cluster_sizes),
        'total_scans': len(df),
        'purity_array': purity_array,
        'size_array': size_array
    }


def main():
    """Main function to run the benchmark"""
    parser = argparse.ArgumentParser(
        description='MS-RT Benchmark for Hyperspec: Calculate N10 and Purity metrics for GPU clustering results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use default parameters
  python msrt_benchmark.py --input data/orbitrap_clusterinfo_gpu_0.265.tsv

  # Specify custom RT window and precursor m/z window
  python msrt_benchmark.py --input data/orbitrap_clusterinfo_gpu_0.265.tsv \\
                           --rt_window 30.0 --precursor_mz_window 0.01

  # Adjust batch size for memory management
  python msrt_benchmark.py --input data/orbitrap_clusterinfo_gpu_0.265.tsv \\
                           --batch_size 5000
        """
    )
    
    parser.add_argument('--input', type=str, required=True,
                       help='Path to GPU cluster info TSV file or parquet file (.parquet or .csv.parquet will be automatically converted)')
    parser.add_argument('--rt_window', type=float, default=30.0,
                       help='Retention time window in seconds (default: 30.0)')
    parser.add_argument('--precursor_mz_window', type=float, default=0.01,
                       help='Precursor m/z window in Da (default: 0.01)')
    parser.add_argument('--batch_size', type=int, default=10000,
                       help='Batch size for multi-threaded processing (affects max memory usage, default: 10000)')
    parser.add_argument('--total_scans', type=int, default=None,
                       help='Total number of scans for N10 calculation (default: number of rows in input file)')
    parser.add_argument('--output_dir', type=str, default='.',
                       help='Output directory for results (default: current directory)')
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)
    
    # Preprocess input file: convert parquet to TSV if needed
    tsv_file = preprocess_input_file(args.input)
    
    # Load results
    # Returns filtered dataframe (for purity) and total scans (for N10, including unclustered)
    gpu_df, total_scans_in_file = load_gpu_results(tsv_file)
    
    # Determine total scans for N10 calculation
    # Note: total scans should include unclustered spectra (cluster == -1)
    # which are filtered out for purity calculation but needed for N10
    if args.total_scans is None:
        # Use total scans from file (includes unclustered spectra)
        total_scans = total_scans_in_file
        print(f"\nUsing total scans from input file (including unclustered): {total_scans:,}")
    else:
        total_scans = args.total_scans
        print(f"\nUsing specified total scans: {total_scans:,}")
    
    # Calculate metrics
    metrics = calculate_metrics(
        gpu_df, 
        total_scans,
        rt_window=args.rt_window,
        precursor_mz_window=args.precursor_mz_window,
        batch_size=args.batch_size
    )
    
    # Print summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print(f"Total Clusters: {metrics['total_clusters']:,}")
    print(f"Total Scans: {metrics['total_scans']:,}")
    print(f"Completeness (N10 metric): {metrics['n10']:,}")
    print(f"Weighted Average Purity: {metrics['weighted_avg_purity']:.6f}")
    print(f"Simple Average Purity: {metrics['simple_avg_purity']:.6f}")
    print(f"{'='*80}\n")
    
    # Save results to file
    output_file = os.path.join(args.output_dir, "msrt_benchmark_results.txt")
    with open(output_file, 'w') as f:
        f.write("MS-RT Benchmark Results\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Input File: {args.input}\n")
        if tsv_file != args.input:
            f.write(f"Converted TSV File: {tsv_file}\n")
        f.write(f"RT Window: {args.rt_window} seconds\n")
        f.write(f"Precursor m/z Window: {args.precursor_mz_window} Da\n")
        f.write(f"Batch Size: {args.batch_size}\n")
        f.write(f"Total Scans: {total_scans:,}\n\n")
        f.write(f"Total Clusters: {metrics['total_clusters']:,}\n")
        f.write(f"Total Scans (in clusters): {metrics['total_scans']:,}\n")
        f.write(f"Completeness (N10 metric): {metrics['n10']:,}\n")
        f.write(f"Weighted Average Purity: {metrics['weighted_avg_purity']:.6f}\n")
        f.write(f"Simple Average Purity: {metrics['simple_avg_purity']:.6f}\n")
    
    print(f"Results saved to: {output_file}")
    print("Done!")


if __name__ == "__main__":
    main()
