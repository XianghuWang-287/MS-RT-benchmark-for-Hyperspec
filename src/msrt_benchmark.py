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


def extract_filename_and_scan_from_identifier(identifier):
    """
    Extract filename and scan number from identifier string.
    
    Handles formats like:
    - "filename.mzML:scan_number" -> (filename.mzML, scan_number)
    - "filename:scan_number" -> (filename, scan_number)
    - "filename_scan_number" -> (filename, scan_number) if scan is numeric
    
    Args:
        identifier: Identifier string that may contain filename and scan
        
    Returns:
        tuple: (filename, scan_number) where scan_number is int or None
    """
    identifier_str = str(identifier)
    
    # Try format: "filename:scan" (most common)
    if ':' in identifier_str:
        parts = identifier_str.rsplit(':', 1)  # Split from right to handle filenames with colons
        if len(parts) == 2:
            filename = parts[0]
            scan_str = parts[1]
            try:
                scan_number = int(scan_str)
                return filename, scan_number
            except ValueError:
                # If scan part is not numeric, return as-is
                return identifier_str, None
    
    # Try format: "filename_scan" where scan is numeric at the end
    # This is a fallback, less reliable
    parts = identifier_str.rsplit('_', 1)
    if len(parts) == 2:
        filename_part, scan_str = parts
        try:
            scan_number = int(scan_str)
            # Only use this if scan looks reasonable (e.g., > 0)
            if scan_number > 0:
                return filename_part, scan_number
        except ValueError:
            pass
    
    # If no pattern matches, return identifier as filename and None for scan
    return identifier_str, None


def load_gpu_results(gpu_file, extract_identifier=False):
    """
    Load GPU clustering results
    
    Args:
        gpu_file: Path to GPU results TSV file
        extract_identifier: If True, extract filename and scan from identifier column
    
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
    
    # Extract filename and scan from identifier if requested
    if extract_identifier:
        print("Extracting filename and scan number from identifier column...")
        print("Note: Original scan column will be ignored, using extracted scan values only.")
        extracted = df['identifier'].apply(extract_filename_and_scan_from_identifier)
        df['extracted_filename'] = extracted.apply(lambda x: x[0])
        df['extracted_scan'] = extracted.apply(lambda x: x[1])
        
        # Use extracted filename and scan (completely replace, ignore original scan column)
        df['filename'] = df['extracted_filename']
        df['scan'] = df['extracted_scan']
        
        # Check if any scans are None or invalid
        if df['scan'].isna().any():
            print(f"Warning: {df['scan'].isna().sum()} identifiers could not be parsed for scan number")
            # Fill missing scans with sequential numbers starting from 1
            missing_mask = df['scan'].isna()
            start_idx = 1
            df.loc[missing_mask, 'scan'] = range(start_idx, start_idx + missing_mask.sum())
        
        df['base_filename'] = df['filename']
        print(f"Extracted filenames: {df['filename'].nunique()} unique files")
        print(f"Scan numbers range: {df['scan'].min()} to {df['scan'].max()}")
    else:
        # GPU identifier format: 000011026_RA3_01_6002 (this is the base filename)
        # scan column is read directly from the file (may be -1)
        # For purity calculation, we'll use identifier as filename
        df['base_filename'] = df['identifier'].astype(str)
        df['filename'] = df['identifier'].astype(str) + '.mzML'  # Add .mzML for consistency
        
        # Use scan column directly from file
        # scan column is read as-is from the TSV file
        # If all scans are -1, assign sequential numbers starting from 1
        if 'scan' in df.columns and (df['scan'] == -1).all():
            print("All scan values are -1, assigning sequential scan numbers starting from 1...")
            df['scan'] = range(1, len(df) + 1)
    
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


def optimized_create_matching_network(cluster, method, rt_window=30.0, precursor_mz_window=0.01, ignore_identifier=False):
    """
    Create matching network for a cluster (optimized version)
    
    Args:
        cluster: DataFrame with cluster data
        method: Method name (e.g., 'falcon')
        rt_window: Retention time window in seconds (default: 30.0)
        precursor_mz_window: Precursor m/z window in Da (default: 0.01)
        ignore_identifier: If True, treat all spectra in cluster as from same file (default: False)
    """
    G = nx.Graph()

    # Precompute the node names and add them to the graph
    # Use row index to ensure each spectrum has a unique node name
    # Even if filename and scan are the same, different spectra (with different m/z or RT) should be separate nodes
    node_attrs = {}
    specs = []
    for idx, (_, row) in enumerate(cluster.iterrows()):
        # Create unique node name using index to handle cases where multiple spectra have same filename_scan
        node_name = f"{row[method_dic[method]['filename']]}_{row[method_dic[method]['scan']]}_{idx}"
        node_attrs[node_name] = {
            "filename": row[method_dic[method]['filename']]
        }
        # Store node_name, filename, mass, rt_time for edge creation
        specs.append((
            node_name,  # node_name for graph (unique per spectrum)
            row[method_dic[method]['filename']],  # filename for comparison (separate to avoid parsing)
            row[method_dic[method]['mass']],
            row[method_dic[method]['rt_time']]
        ))
    G.add_nodes_from(node_attrs.items())

    # Create edges based on conditions
    # RT tolerance: configurable (default: 30 seconds)
    # Precursor m/z tolerance: configurable (default: 0.01 Da)
    # 
    # PERFORMANCE OPTIMIZATION: Group by filename first to avoid generating unnecessary combinations
    # For clusters with multiple files, this dramatically reduces the number of combinations to check
    # Example: 10k spectra in 10 files -> only check 10 × C(1k,2) instead of C(10k,2)
    # 
    # Add a global/procedural variable to control tqdm/progress bar usage
    # By default, progress bar is off, unless set elsewhere by the user
    USE_PROGRESS_BAR = globals().get("USE_PROGRESS_BAR", False)

    if ignore_identifier:
        # When ignoring identifier, all spectra are treated as from same file
        # PERFORMANCE OPTIMIZATION: Use sliding window to avoid checking all combinations
        # Sort by m/z first, then only check nearby spectra within the m/z window
        # This dramatically reduces combinations for large clusters
        # Example: 10k spectra -> from 50M to ~100k combinations (500x speedup)
        
        # Sort by m/z (index 2) for efficient sliding window
        specs_sorted = sorted(specs, key=lambda x: x[2])  # x[2] is m/z
        
        # Use sliding window technique: for each spectra, only check nearby ones
        def generate_edges_sliding_window():
            n = len(specs_sorted)
            for i in range(n):
                spec1 = specs_sorted[i]
                mz1 = spec1[2]
                rt1 = spec1[3]
                
                # Find the range of spectra with m/z within window
                # Since sorted, we can use early termination
                j = i + 1
                while j < n:
                    spec2 = specs_sorted[j]
                    mz2 = spec2[2]
                    
                    # Early termination: if m/z difference exceeds window, stop checking
                    if mz2 - mz1 > precursor_mz_window:
                        break
                    
                    # Check both m/z and RT windows
                    if abs(mz1 - mz2) <= precursor_mz_window:
                        rt2 = spec2[3]
                        if abs(rt1 - rt2) <= rt_window:
                            yield (spec1[0], spec2[0])
                    
                    j += 1
        
        # Generate edges using sliding window
        edges = generate_edges_sliding_window()
        
        if USE_PROGRESS_BAR:
            from tqdm import tqdm as _tqdm
            # Estimate total combinations for progress bar
            # Actual will be much less due to sliding window optimization
            n = len(specs)
            total_combinations = n * (n - 1) // 2
            edges = _tqdm(edges, total=total_combinations, desc="Building cluster network")
    else:
        # OPTIMIZATION: Group by filename first, then only generate combinations within each group
        # This avoids generating combinations between different files
        from collections import defaultdict
        
        # Group specs by filename
        specs_by_file = defaultdict(list)
        for spec in specs:
            filename = spec[1]  # spec[1] is the filename
            specs_by_file[filename].append(spec)
        
        # Generate edges only within each file group
        total_combinations = 0
        edge_generators = []
        
        for filename, file_specs in specs_by_file.items():
            if len(file_specs) < 2:
                continue  # Skip files with less than 2 spectra
            
            # Generate combinations only within this file
            file_comb_iter = combinations(file_specs, 2)
            total_combinations += len(file_specs) * (len(file_specs) - 1) // 2
            
            # Create edge generator for this file
            file_edges = (
                (spec1[0], spec2[0]) for spec1, spec2 in file_comb_iter
                if abs(spec1[2] - spec2[2]) <= precursor_mz_window 
                   and abs(spec1[3] - spec2[3]) <= rt_window
            )
            edge_generators.append(file_edges)
        
        # Combine all edge generators
        from itertools import chain
        edges = chain(*edge_generators)
        
        if USE_PROGRESS_BAR:
            from tqdm import tqdm as _tqdm
            edges = _tqdm(edges, total=total_combinations, desc="Building cluster network")
    
    # add_edges_from can accept a generator, which is memory efficient
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


def process_cluster(cluster_id, cluster_data, method='falcon', rt_window=30.0, precursor_mz_window=0.01, ignore_identifier=False):
    """
    Process a single cluster to calculate purity (MS-RT method)
    
    Args:
        cluster_id: Cluster ID
        cluster_data: DataFrame with cluster data
        method: Method name (e.g., 'falcon')
        rt_window: Retention time window in seconds
        precursor_mz_window: Precursor m/z window in Da
        ignore_identifier: If True, treat all spectra in cluster as from same file (default: False)
    """
    if len(cluster_data) == 1:
        return (1.0, 1)  # Singleton cluster has purity 1
    
    G = optimized_create_matching_network(cluster_data, method, rt_window, precursor_mz_window, ignore_identifier)
    
    # Calculate the count of each filename in the cluster
    # If ignore_identifier is True, treat all as from same file
    if ignore_identifier:
        # All spectra are treated as from the same file
        file_counts = Counter({'same_file': len(cluster_data)})
        # Find the largest connected component in the entire graph
        # Use generator instead of converting to list to save memory for large graphs
        components = nx.connected_components(G)
        largest_component_size = 1
        for comp in components:
            comp_size = len(comp)
            if comp_size > largest_component_size:
                largest_component_size = comp_size
        max_component_sizes = {'same_file': largest_component_size}
    else:
        file_counts = Counter(cluster_data[method_dic[method]['filename']])
        max_component_sizes = calculate_max_component_per_file(G)
    
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


def falcon_purity(cluster_results, batch_size=10000, rt_window=30.0, precursor_mz_window=0.01, ignore_identifier=False):
    """
    Calculate purity for each cluster using network-based approach.
    Uses batching to avoid loading all clusters into memory at once.
    
    Args:
        cluster_results: DataFrame with cluster assignments
        batch_size: Number of clusters to process in each batch (default: 10000)
        rt_window: Retention time window in seconds (default: 30.0)
        precursor_mz_window: Precursor m/z window in Da (default: 0.01)
        ignore_identifier: If True, treat all spectra in each cluster as from same file (default: False)
    
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
    if ignore_identifier:
        print("Ignoring identifier: treating all spectra in each cluster as from same file")
    
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
                                     rt_window=rt_window, precursor_mz_window=precursor_mz_window,
                                     ignore_identifier=ignore_identifier)
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
                                 rt_window=rt_window, precursor_mz_window=precursor_mz_window,
                                 ignore_identifier=ignore_identifier)
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


def calculate_metrics(df, total_scans, rt_window=30.0, precursor_mz_window=0.01, batch_size=10000, ignore_identifier=False):
    """
    Calculate N10 and Purity metrics for a clustering result
    
    Args:
        df: DataFrame with cluster assignments
        total_scans: Total number of scans for N10 calculation
        rt_window: Retention time window in seconds (default: 30.0)
        precursor_mz_window: Precursor m/z window in Da (default: 0.01)
        batch_size: Batch size for processing (default: 10000)
        ignore_identifier: If True, treat all spectra in each cluster as from same file (default: False)
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
                                          precursor_mz_window=precursor_mz_window,
                                          ignore_identifier=ignore_identifier)
    
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
    parser.add_argument('--rt_window', type=float, default=None,
                       help='Retention time window value; unit is controlled by --rt_unit '
                            '(default: 30.0 seconds or 0.5 minutes)')
    parser.add_argument('--rt_unit', type=str, choices=['seconds', 'minutes'], default='seconds',
                       help='Unit for retention time window: "seconds" or "minutes" (default: seconds)')
    parser.add_argument('--precursor_mz_window', type=float, default=0.01,
                       help='Precursor m/z window in Da (default: 0.01)')
    parser.add_argument('--batch_size', type=int, default=10000,
                       help='Batch size for multi-threaded processing (affects max memory usage, default: 10000)')
    parser.add_argument('--total_scans', type=int, default=None,
                       help='Total number of scans for N10 calculation (default: number of rows in input file)')
    parser.add_argument('--ignore_identifier', action='store_true',
                       help='Ignore identifier: treat all spectra in each cluster as from same file (default: False)')
    parser.add_argument('--extract_identifier', action='store_true',
                       help='Extract filename and scan number from identifier column (e.g., "file.mzML:123" -> filename="file.mzML", scan=123)')
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
    gpu_df, total_scans_in_file = load_gpu_results(tsv_file, extract_identifier=args.extract_identifier)
    
    # Convert retention_time in data if unit is minutes
    # If rt_unit is minutes, data retention_time values are in minutes, convert to seconds
    if args.rt_unit == 'minutes':
        print(f"Converting retention_time from minutes to seconds (multiplying by 60)...")
        gpu_df['retention_time'] = gpu_df['retention_time'] * 60.0
    
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
    
    # Determine RT window in seconds
    # - If unit is seconds and no value is provided: default to 30.0 seconds
    # - If unit is minutes and no value is provided: default to 0.5 minutes (30.0 seconds)
    # Note: After conversion, all retention_time values in data are in seconds
    if args.rt_unit == 'minutes':
        if args.rt_window is None:
            rt_window_seconds = 0.5 * 60.0
        else:
            rt_window_seconds = args.rt_window * 60.0
    else:
        if args.rt_window is None:
            rt_window_seconds = 30.0
        else:
            rt_window_seconds = args.rt_window

    rt_window_minutes = rt_window_seconds / 60.0
    
    # Calculate metrics
    metrics = calculate_metrics(
        gpu_df, 
        total_scans,
        rt_window=rt_window_seconds,
        precursor_mz_window=args.precursor_mz_window,
        batch_size=args.batch_size,
        ignore_identifier=args.ignore_identifier
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
        f.write(f"RT Window: {rt_window_seconds} seconds ({rt_window_minutes:.3f} minutes)\n")
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
