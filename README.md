# MS-RT Benchmark for Hyperspec

A Python tool for evaluating GPU clustering results using MS-RT (Mass Spectrometry Retention Time) validation metrics^[1]. This tool calculates N10 (completeness proxy) and Purity metrics to assess the quality of clustering results from GPU-based clustering algorithms. The MS-RT method provides a practical approach for evaluating MS/MS clustering performance in metabolomics, where false discovery rate (FDR)-controlled identifications are not currently practical.

## Overview

This benchmark tool processes GPU clustering results in TSV format and calculates two key metrics:

1. **Completeness (N10)**: The cluster size at which 10% of total scans are covered. This metric helps assess the distribution of cluster sizes and identifies the threshold for large clusters.

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
pyarrow  # Required for parquet file support
```

Install dependencies using pip:

```bash
pip install pandas numpy networkx tqdm pyarrow
```

Or install from requirements.txt:

```bash
pip install -r requirements.txt
```

### Input File Format

The tool supports two input formats:

1. **TSV files** (`.tsv`): Direct input format
2. **Parquet files** (`.parquet` or `.csv.parquet`): Will be automatically converted to TSV during preprocessing

The input file (TSV or parquet) should contain the following columns:

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

**Note on Parquet Files**: If you provide a `.parquet` or `.csv.parquet` file, the tool will automatically:
1. Detect the file format
2. Convert it to TSV format using the `convert_to_tsv.py` script (or direct pandas conversion)
3. Save the converted TSV file in the same directory with the naming pattern: `<dataset>_clusterinfo_gpu.tsv`
4. Use the converted TSV file for benchmarking

The converted TSV file will be reused if it already exists, so you don't need to reconvert the same file multiple times.

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
# Using TSV file
python src/msrt_benchmark.py --input data/orbitrap_clusterinfo_gpu.tsv

# Using parquet file (will be automatically converted)
python src/msrt_benchmark.py --input data/orbitrap_clustering_results_1k.csv.parquet
```

This will use all default parameters:
- RT window: 30.0 seconds
- Precursor m/z window: 0.01 Da
- Batch size: 10000

### Input File Formats

#### Using TSV Files

If you already have TSV files, you can use them directly:

```bash
python src/msrt_benchmark.py --input data/orbitrap_clusterinfo_gpu.tsv
```

#### Using Parquet Files

The tool automatically detects and converts parquet files:

```bash
# Single parquet file
python src/msrt_benchmark.py --input data/orbitrap_clustering_results_1k.csv.parquet

# The tool will automatically:
# 1. Detect it's a parquet file
# 2. Convert to TSV (saved as orbitrap_clusterinfo_gpu.tsv)
# 3. Run the benchmark on the converted file
```

**Batch Conversion**: If you have multiple parquet files to convert, you can use the `convert_to_tsv.py` script directly for batch processing:

```bash
# Convert a single file
python src/convert_to_tsv.py data/orbitrap_clustering_results_1k.csv.parquet data/orbitrap_clusterinfo_gpu.tsv

# Convert with auto-generated output name
python src/convert_to_tsv.py data/orbitrap_clustering_results_1k.csv.parquet

# Then use the converted TSV files
python src/msrt_benchmark.py --input data/orbitrap_clusterinfo_gpu.tsv
```

The `convert_to_tsv.py` script can be useful for:
- Batch converting multiple parquet files before benchmarking
- Pre-processing files to avoid conversion overhead during benchmarking
- Converting files for other purposes

### Advanced Usage

#### Custom Retention Time Window

Adjust the retention time window for matching spectra within clusters:

```bash
python src/msrt_benchmark.py --input data/orbitrap_clusterinfo_gpu.tsv \
                             --rt_window 30.0
```

The RT window is specified in seconds. A larger window allows spectra with more distant retention times to be considered as matches, which may increase connectivity in the matching network but could also reduce purity scores.

#### Custom Precursor m/z Window

Adjust the precursor m/z tolerance for matching spectra:

```bash
python src/msrt_benchmark.py --input data/orbitrap_clusterinfo_gpu.tsv \
                             --precursor_mz_window 0.01
```

The precursor m/z window is specified in Daltons (Da). A larger window allows spectra with more different precursor m/z values to be considered as matches.

#### Adjust Batch Size for Memory Management

Control memory usage by adjusting the batch size:

```bash
python src/msrt_benchmark.py --input data/orbitrap_clusterinfo_gpu.tsv \
                             --batch_size 5000
```

**Important Notes on Batch Size:**
- Smaller batch sizes use less memory but may be slightly slower
- Larger batch sizes use more memory but may be faster
- Default value (10000) is optimized for most systems
- If you encounter memory errors, reduce the batch size (e.g., 5000 or 2000)
- If you have abundant memory, increasing batch size (e.g., 20000) may improve performance

#### Specify Total Scans for N10 Calculation

By default, the total number of scans is automatically inferred from the input file (the number of rows in the file after filtering out noise clusters with `cluster == -1`). The tool will display: `Using total scans from input file: X`.

**When to specify `--total_scans`:**
- If you want to include noise clusters (with `cluster == -1`) that were filtered out in the N10 calculation
- If you have a known total number of scans from the original dataset before clustering
- If you want to compare N10 values across different clustering runs using the same total scan count

**Example:**
```bash
# Default: uses number of rows in input file (1,237,710 in this example)
python src/msrt_benchmark.py --input data/orbitrap_clusterinfo_gpu.tsv

# Specify a different total (e.g., including noise clusters)
python src/msrt_benchmark.py --input data/orbitrap_clusterinfo_gpu.tsv \
                             --total_scans 1500000
```

**Note:** The N10 metric is calculated based on the cumulative coverage of scans. Specifying a larger `--total_scans` value will result in a higher N10 value, as it requires more scans to reach the 10% coverage threshold.

#### Custom Output Directory

Specify where to save the results:

```bash
python src/msrt_benchmark.py --input data/orbitrap_clusterinfo_gpu.tsv \
                             --output_dir results/
```

### Complete Example

A complete example with all parameters specified:

```bash
python src/msrt_benchmark.py \
    --input data/orbitrap_clusterinfo_gpu.tsv \
    --rt_window 30.0 \
    --precursor_mz_window 0.01 \
    --batch_size 10000 \
    --total_scans 1237710 \
    --output_dir results/
```

## Command Line Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--input` | string | **Required** | Path to GPU cluster info TSV file or parquet file (.parquet or .csv.parquet will be automatically converted) |
| `--rt_window` | float | 30.0 | Retention time window in seconds |
| `--precursor_mz_window` | float | 0.01 | Precursor m/z window in Da |
| `--batch_size` | int | 10000 | Batch size for multi-threaded processing (affects max memory usage) |
| `--total_scans` | int | None | Total number of scans for N10 calculation (default: automatically inferred from input file, displays "Using total scans from input file: X") |
| `--output_dir` | string | . | Output directory for results |

## Output

The tool generates two types of output:

### 1. Console Output

The tool prints progress information and a summary to the console:

```
Loading GPU results from ./data/orbitrap_clusterinfo_gpu.tsv...
Loaded 1237710 GPU scan assignments
GPU clusters: 904976

Using total scans from input file: 1,237,710

================================================================================
Calculating metrics for GPU clustering results
================================================================================
Calculating N10...
N10 value: 22
Calculating purity (this may take a while)...
Processing 904,976 clusters in batches of 10,000...
Using RT window: 30.0 seconds, Precursor m/z window: 0.01 Da
Calculating purity: 100%|██████████████████████████████████████████████████████████████████████████████████████████████████████| 904976/904976 [03:49<00:00, 3937.05it/s]
Processed 904,976 clusters successfully.
Weighted Average Purity: 0.979450
Simple Average Purity: 0.995611

================================================================================
SUMMARY
================================================================================
Total Clusters: 904,976
Total Scans: 1,237,710
Completeness (N10 metric): 22
Weighted Average Purity: 0.979450
Simple Average Purity: 0.995611
================================================================================

Results saved to: ./msrt_benchmark_results.txt
Done!
```

### 2. Results File

A text file (`msrt_benchmark_results.txt`) is saved to the output directory containing:

- Input file path (and converted TSV file path if parquet was used)
- Parameter settings (RT window, precursor m/z window, batch size)
- Total scans (used for N10 calculation)
- Total clusters
- Total scans in clusters (after filtering noise)
- Completeness (N10 metric)
- Weighted average purity
- Simple average purity

**Example results file content:**
```
MS-RT Benchmark Results
================================================================================

Input File: data/orbitrap_clusterinfo_gpu.tsv
RT Window: 30.0 seconds
Precursor m/z Window: 0.01 Da
Batch Size: 10000
Total Scans: 1237710

Total Clusters: 904,976
Total Scans (in clusters): 1,237,710
Completeness (N10 metric): 22
Weighted Average Purity: 0.979450
Simple Average Purity: 0.995611
```

## Understanding the Metrics

### N10 (Completeness Metric)

N10 represents the cluster size threshold at which 10% of all scans are covered. This metric serves as a proxy for clustering completeness, which evaluates the extent to which all members of a true class are assigned to the same cluster. In metabolomics, where false discovery rate (FDR)-controlled identifications are not currently practical, N10 provides a practical alternative to traditional completeness metrics used in proteomics.

**Relationship to Completeness:**
- **Completeness** measures an algorithm's ability to ensure that all elements belonging to a particular group are clustered together, reflecting the clustering method's capacity to capture the entirety of natural groupings within the data
- **N10** serves as a proxy for completeness by identifying the cluster size threshold that covers a significant portion (10%) of the data, indicating how well the clustering algorithm consolidates related spectra into large, comprehensive clusters
- Higher N10 values suggest better completeness, as they indicate that a substantial portion of scans are successfully grouped into large clusters, demonstrating the algorithm's effectiveness in capturing complete groupings

**Interpretation:**
- **Lower N10**: Indicates that many scans are in relatively small clusters, suggesting lower completeness
- **Higher N10**: Indicates that a significant portion of scans are in large clusters, suggesting higher completeness and better consolidation of related spectra

**Typical N10 Values:**
- For datasets with ~1M scans: N10 values typically range from 10-100, depending on clustering parameters
- For datasets with ~10M scans: With appropriate clustering parameters (e.g., eps ~0.25 for Falcon), N10 values in the thousands (e.g., 1,000-5,000) are expected
- **N10 values that are too small** (e.g., < 100 for 10M scans) may indicate that the clustering algorithm is splitting related spectra into too many small clusters, suggesting the need to adjust clustering parameters (e.g., increase eps for Falcon) to achieve better consolidation

This metric has been validated against traditional database search-based completeness metrics in proteomics datasets, demonstrating consistent relative performance ordering of clustering tools^[1].

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

**Typical Purity Values:**
- **Typical purity range**: Well-performing clustering algorithms typically achieve purity values in the range of **0.97-0.99** (weighted average purity)
- **Purity < 0.97**: May indicate clustering quality issues, such as:
  - Over-clustering: Related spectra are being split into separate clusters
  - Under-clustering: Unrelated spectra are being grouped together
  - Parameter misconfiguration: RT window or precursor m/z window may need adjustment
- **Purity > 0.99**: Excellent clustering quality, indicating that most clusters contain well-connected, related spectra

**Note:** The actual purity value depends on the clustering algorithm, dataset characteristics, and parameter settings. Values consistently below 0.97 may warrant investigation of clustering parameters or algorithm selection.

## Performance Considerations

### Memory Usage

The tool uses batch processing to manage memory efficiently. The `--batch_size` parameter directly affects maximum memory usage:

- **Default (10000)**: Suitable for most systems with 8-16 GB RAM
- **Smaller (2000-5000)**: Use if you encounter memory errors
- **Larger (20000+)**: Use if you have abundant memory (32+ GB RAM)

### Processing Time

Processing time depends on:
- Number of clusters (more clusters = longer processing)
- Cluster sizes (larger clusters take more time to process due to more network edge calculations)
- Batch size (smaller batches may be slightly slower due to overhead)
- Number of CPU cores (parallel processing is used)
- RT window and precursor m/z window (larger windows create more edges, increasing computation)

**Performance Benchmarks:**
Based on actual test runs with default parameters (RT window: 30.0 seconds, precursor m/z window: 0.01 Da, batch size: 10,000):
- **~900K clusters, ~1.2M scans**: ~4 minutes (processing rate: ~3,900 clusters/second)
- Processing rate typically ranges from 2,000-5,000 clusters/second depending on cluster size distribution and hardware

**Typical processing times (estimated):**
- Small datasets (< 100K clusters): 10-30 seconds
- Medium datasets (100K - 500K clusters): 1-5 minutes
- Medium-large datasets (500K - 1M clusters): 3-10 minutes
- Large datasets (1M - 5M clusters): 10-60 minutes
- Very large datasets (> 5M clusters): 1-4 hours

**Note:** Actual processing time can vary significantly based on:
- Average cluster size (larger clusters require more network computations)
- Hardware specifications (CPU cores, memory speed)
- System load during processing

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
--input /full/path/to/data/orbitrap_clusterinfo_gpu.tsv
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
├── requirements.txt       # Python dependencies
├── .gitignore            # Git ignore file
├── src/
│   ├── msrt_benchmark.py  # Main benchmark script
│   └── convert_to_tsv.py  # Utility script for converting parquet to TSV
└── data/                  # Place your input TSV or parquet files here
```

## Additional Tools

### convert_to_tsv.py

A utility script for converting parquet files to TSV format. This script can be used independently for batch conversion or is automatically called by `msrt_benchmark.py` when parquet files are detected.

**Usage:**

```bash
# Convert with auto-generated output filename
python src/convert_to_tsv.py data/orbitrap_clustering_results_1k.csv.parquet

# Specify output filename
python src/convert_to_tsv.py data/orbitrap_clustering_results_1k.csv.parquet data/orbitrap_clusterinfo_gpu.tsv

# Convert with full path
python src/convert_to_tsv.py /path/to/orbitrap_clustering_results_1k.csv.parquet
```

**Features:**
- Automatically extracts dataset name from file path (e.g., MSV numbers)
- Generates output filename: `<dataset>_clusterinfo_gpu.tsv`
- Shows file size and column information after conversion
- Can be used for batch processing multiple parquet files

## License

[Specify your license here]

## Citation

If you use this tool in your research, please cite the MS-RT method paper:

**Wang, X.**, El Abiead, Y., Acharya, D. D., Brown, C. J., Clevenger, K., Hu, J., Kretsch, A., Menegatti, C., Xiong, Q., Bittremieux, W., & Wang, M. (2025). MS-RT: A Method for Evaluating MS/MS Clustering Performance for Metabolomics Data. *Journal of Proteome Research*, 24(4), 1778-1790. doi: [10.1021/acs.jproteome.4c00881](https://doi.org/10.1021/acs.jproteome.4c00881)

**BibTeX:**
```bibtex
@article{wang2025msrt,
  title={MS-RT: A Method for Evaluating MS/MS Clustering Performance for Metabolomics Data},
  author={Wang, Xianghu and El Abiead, Yasin and Acharya, Deepa D and Brown, Christopher J and Clevenger, Ken and Hu, Jie and Kretsch, Ashley and Menegatti, Carla and Xiong, Quanbo and Bittremieux, Wout and Wang, Mingxun},
  journal={Journal of Proteome Research},
  volume={24},
  number={4},
  pages={1778--1790},
  year={2025},
  publisher={ACS Publications},
  doi={10.1021/acs.jproteome.4c00881}
}
```

**Reference:**
[1] Wang, X., et al. (2025). MS-RT: A Method for Evaluating MS/MS Clustering Performance for Metabolomics Data. *Journal of Proteome Research*, 24(4), 1778-1790.

## Contact

[Add contact information]

## Version History

- **v1.0.0**: Initial release
  - Support for GPU clustering results
  - Configurable RT window, precursor m/z window, and batch size
  - N10 and Purity metric calculation
