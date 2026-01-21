#!/usr/bin/env python3
"""
Convert clustering results from parquet to TSV format
Usage:
    python convert_to_tsv.py [input_file] [output_file]
    
Requirements:
    pip install pandas pyarrow
"""

import sys
import argparse
from pathlib import Path
import re

# Check dependencies
try:
    import pandas as pd
except ImportError:
    print("Error: pandas not found. Install with: pip install pandas pyarrow", file=sys.stderr)
    sys.exit(1)

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

def convert_to_tsv(input_file, output_file=None):
    """Convert parquet file to TSV format"""
    input_path = Path(input_file)
    
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)
    
    # Determine output filename
    if output_file:
        output_path = Path(output_file)
    else:
        # Auto-generate output filename
        dataset_name = extract_dataset_name(input_path)
        output_dir = input_path.parent
        output_path = output_dir / f"{dataset_name}_clusterinfo_gpu.tsv"
    
    print(f"Reading parquet file: {input_path}")
    
    try:
        # Read parquet file
        df = pd.read_parquet(input_path)
        print(f"Loaded {len(df):,} rows, {len(df.columns)} columns")
        
        # Write to TSV
        print(f"Writing TSV file: {output_path}")
        df.to_csv(output_path, sep='\t', index=False)
        
        # Get file size
        file_size = output_path.stat().st_size
        file_size_mb = file_size / (1024 * 1024)
        
        print(f"✓ Conversion completed!")
        print(f"  Output file: {output_path}")
        print(f"  File size: {file_size_mb:.2f} MB")
        print(f"  Rows: {len(df):,}")
        print(f"  Columns: {len(df.columns)}")
        
        # Show column names
        print(f"\nColumns: {', '.join(df.columns.tolist())}")
        
        return output_path
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='Convert clustering results from parquet to TSV format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-generate output filename (MSV000081981_clusterinfo_gpu.tsv)
  python convert_to_tsv.py results.parquet
  
  # Specify output filename
  python convert_to_tsv.py results.parquet output.tsv
  
  # Convert with full path
  python convert_to_tsv.py /path/to/results.parquet
        """
    )
    
    parser.add_argument('input_file', 
                       nargs='?',
                       default='/data/nas-gpu/wang/xianghu/metabolomics_clustering_data/MSV000081981_clustering_results_1k.csv.parquet',
                       help='Input parquet file (default: MSV000081981_clustering_results_1k.csv.parquet)')
    parser.add_argument('output_file', nargs='?', default=None,
                       help='Output TSV file (default: auto-generated as <dataset>_clusterinfo_gpu.tsv)')
    
    args = parser.parse_args()
    
    convert_to_tsv(args.input_file, args.output_file)

if __name__ == '__main__':
    main()
