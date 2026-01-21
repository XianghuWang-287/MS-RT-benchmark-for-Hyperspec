# MS-RT Benchmark for Hyperspec

A Python tool for evaluating GPU clustering results using MS-RT (Mass Spectrometry Retention Time) validation metrics. This tool calculates N10 and Purity metrics to assess the quality of clustering results from GPU-based clustering algorithms.

## Overview

This benchmark tool processes GPU clustering results in TSV format and calculates two key metrics:

1. **N10**: The cluster size at which 10% of total scans are covered. This metric helps assess the distribution of cluster sizes and identifies the threshold for large clusters.

2. **Purity (MS-RT validation)**: A network-based purity metric that evaluates cluster quality by examining the connectivity of spectra within clusters based on retention time and precursor m/z windows. The purity score ranges from 0 to 1, where 1 indicates perfect clustering.

## Features

- **Configurable Parameters**: 
  - Retention time window (default: 30 seconds)
  - Precursor m/z window (default: 0.01 Da)
  - Multi-threaded batch size (default: 10000, affects maximum memory usage)
  
- **Memory Efficient**: Uses batch processing to handle large datasets without excessive memory consumption

- **GPU Results Focused**: Specifically designed for MS-RT benchmark GPU results format

## Requirements

### Python Dependencies

The tool requires the following Python packages:

```
pandas
numpy
networkx
tqdm
```

Install dependencies using pip:

```bash
pip install pandas numpy networkx tqdm
```

### Input File Format

The input TSV file should contain the following columns:

- `identifier`: Sample identifier (e.g., "000011026_RA3_01_6002")
- `precursor_mz`: Precursor m/z value
- `retention_time`: Retention time in seconds
- `cluster`: Cluster assignment (clusters with value -1 are treated as noise and excluded)

Example input file structure:

```
identifier	precursor_mz	retention_time	cluster
000011026_RA3_01_6002	456.7890	1234.56	1
000011026_RA3_01_6002	456.7891	1234.58	1
000011026_RA3_01_6002	789.1234	2345.67	2
```

## Installation

1. Clone or download this repository
2. Ensure all dependencies are installed (see Requirements section)
3. Make the script executable (optional):

```bash
chmod +x src/msrt_benchmark.py
```

## Usage

### Basic Usage

The simplest way to run the benchmark is to specify only the input file:

```bash
python src/msrt_benchmark.py --input data/orbitrap_clusterinfo_gpu_0.265.tsv
```

This will use all default parameters:
- RT window: 30.0 seconds
- Precursor m/z window: 0.01 Da
- Batch size: 10000

### Advanced Usage

#### Custom Retention Time Window

Adjust the retention time window for matching spectra within clusters:

```bash
python src/msrt_benchmark.py --input data/orbitrap_clusterinfo_gpu_0.265.tsv \
                             --rt_window 30.0
```

The RT window is specified in seconds. A larger window allows spectra with more distant retention times to be considered as matches, which may increase connectivity in the matching network but could also reduce purity scores.

#### Custom Precursor m/z Window

Adjust the precursor m/z tolerance for matching spectra:

```bash
python src/msrt_benchmark.py --input data/orbitrap_clusterinfo_gpu_0.265.tsv \
                             --precursor_mz_window 0.01
```

The precursor m/z window is specified in Daltons (Da). A larger window allows spectra with more different precursor m/z values to be considered as matches.

#### Adjust Batch Size for Memory Management

Control memory usage by adjusting the batch size:

```bash
python src/msrt_benchmark.py --input data/orbitrap_clusterinfo_gpu_0.265.tsv \
                             --batch_size 5000
```

**Important Notes on Batch Size:**
- Smaller batch sizes use less memory but may be slightly slower
- Larger batch sizes use more memory but may be faster
- Default value (10000) is optimized for most systems
- If you encounter memory errors, reduce the batch size (e.g., 5000 or 2000)
- If you have abundant memory, increasing batch size (e.g., 20000) may improve performance

#### Specify Total Scans for N10 Calculation

By default, the total number of scans is inferred from the input file. If you need to use a different value (e.g., including noise clusters that were filtered out):

```bash
python src/msrt_benchmark.py --input data/orbitrap_clusterinfo_gpu_0.265.tsv \
                             --total_scans 7342436
```

#### Custom Output Directory

Specify where to save the results:

```bash
python src/msrt_benchmark.py --input data/orbitrap_clusterinfo_gpu_0.265.tsv \
                             --output_dir results/
```

### Complete Example

A complete example with all parameters specified:

```bash
python src/msrt_benchmark.py \
    --input data/orbitrap_clusterinfo_gpu_0.265.tsv \
    --rt_window 30.0 \
    --precursor_mz_window 0.01 \
    --batch_size 10000 \
    --total_scans 7342436 \
    --output_dir results/
```

## Command Line Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--input` | string | **Required** | Path to GPU cluster info TSV file |
| `--rt_window` | float | 30.0 | Retention time window in seconds |
| `--precursor_mz_window` | float | 0.01 | Precursor m/z window in Da |
| `--batch_size` | int | 10000 | Batch size for multi-threaded processing (affects max memory usage) |
| `--total_scans` | int | None | Total number of scans for N10 calculation (default: number of rows in input file) |
| `--output_dir` | string | . | Output directory for results |

## Output

The tool generates two types of output:

### 1. Console Output

The tool prints progress information and a summary to the console:

```
Loading GPU results from data/orbitrap_clusterinfo_gpu_0.265.tsv...
Loaded 7342436 GPU scan assignments
GPU clusters: 1234567

================================================================================
Calculating metrics for GPU clustering results
================================================================================
Calculating N10...
N10 value: 12345
Calculating purity (this may take a while)...
Processing 1,234,567 clusters in batches of 10,000...
Using RT window: 30.0 seconds, Precursor m/z window: 0.01 Da
Calculating purity: 100%|████████████████████| 1234567/1234567 [10:23<00:00, 1987.45it/s]
Processed 1,234,567 clusters successfully.
Weighted Average Purity: 0.923456
Simple Average Purity: 0.912345

================================================================================
SUMMARY
================================================================================
Total Clusters: 1,234,567
Total Scans: 7,342,436
N10: 12,345
Weighted Average Purity: 0.923456
Simple Average Purity: 0.912345
================================================================================

Results saved to: ./msrt_benchmark_results.txt
Done!
```

### 2. Results File

A text file (`msrt_benchmark_results.txt`) is saved to the output directory containing:

- Input file path
- Parameter settings (RT window, precursor m/z window, batch size)
- Total scans
- Total clusters
- N10 value
- Weighted average purity
- Simple average purity

## Understanding the Metrics

### N10

N10 represents the cluster size threshold at which 10% of all scans are covered. This metric helps understand the distribution of cluster sizes:

- **Lower N10**: Indicates that many scans are in relatively small clusters
- **Higher N10**: Indicates that a significant portion of scans are in large clusters

### Purity (MS-RT Validation)

Purity is calculated using a network-based approach:

1. For each cluster, a matching network is constructed where nodes represent spectra
2. Edges connect spectra that:
   - Come from the same file (identifier)
   - Have precursor m/z values within the specified window
   - Have retention times within the specified window
3. The purity for each file in the cluster is the fraction of its spectra in the largest connected component
4. The overall cluster purity is a weighted average across all files in the cluster

**Interpretation:**
- **Purity = 1.0**: All spectra in the cluster form a single connected component (perfect clustering)
- **Purity < 1.0**: Some spectra are disconnected, suggesting potential clustering errors
- **Higher purity**: Generally indicates better clustering quality

## Performance Considerations

### Memory Usage

The tool uses batch processing to manage memory efficiently. The `--batch_size` parameter directly affects maximum memory usage:

- **Default (10000)**: Suitable for most systems with 8-16 GB RAM
- **Smaller (2000-5000)**: Use if you encounter memory errors
- **Larger (20000+)**: Use if you have abundant memory (32+ GB RAM)

### Processing Time

Processing time depends on:
- Number of clusters (more clusters = longer processing)
- Cluster sizes (larger clusters take more time to process)
- Batch size (smaller batches may be slightly slower due to overhead)
- Number of CPU cores (parallel processing is used)

Typical processing times:
- Small datasets (< 100K clusters): Seconds to minutes
- Medium datasets (100K - 1M clusters): Minutes to tens of minutes
- Large datasets (> 1M clusters): Tens of minutes to hours

## Troubleshooting

### Memory Errors

If you encounter memory errors (e.g., "MemoryError" or "Killed"):

1. Reduce the batch size:
   ```bash
   --batch_size 2000
   ```

2. Ensure no other memory-intensive processes are running

3. Consider processing on a machine with more RAM

### File Not Found Errors

Ensure the input file path is correct. Use absolute paths if relative paths don't work:

```bash
--input /full/path/to/data/orbitrap_clusterinfo_gpu_0.265.tsv
```

### Slow Processing

If processing is very slow:

1. Check that you have multiple CPU cores available (the tool uses multiprocessing)
2. Consider increasing batch size if you have available memory
3. Verify the input file format is correct

## File Structure

```
MS-RT benchmark for Hyperspec/
├── README.md              # This file
├── src/
│   └── msrt_benchmark.py  # Main benchmark script
└── data/                  # Place your input TSV files here
```

## License

[Specify your license here]

## Citation

If you use this tool in your research, please cite:

[Add citation information]

## Contact

[Add contact information]

## Version History

- **v1.0.0**: Initial release
  - Support for GPU clustering results
  - Configurable RT window, precursor m/z window, and batch size
  - N10 and Purity metric calculation
