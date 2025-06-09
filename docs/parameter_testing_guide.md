# DNA Sequence Reconstruction - Parameter Testing Guide

A comprehensive guide for testing different parameters with the quality analysis script to optimize DNA sequence reconstruction performance.

## Table of Contents

- [Quick Start](#quick-start)
- [Command Structure](#command-structure)
- [Analysis Types](#analysis-types)
- [Parameter Testing Examples](#parameter-testing-examples)
- [Parameter Reference](#parameter-reference)
- [Best Practices](#best-practices)
- [Reading Results](#reading-results)
- [Troubleshooting](#troubleshooting)

---

## Quick Start

```bash
# Change to scripts directory
cd scripts

# Run a quick test to verify everything works
python quality_analysis.py --analysis error --k 10 --n 50 --length 100 --error-rates 0.0 0.1 --repeats 1 --output ../data/results/quick_test

# Check results
cat ../data/results/quick_test/quality_analysis_report.txt
```

---

## Command Structure

### Basic Template

```bash
python quality_analysis.py --analysis [TYPE] [PARAMETERS] --output ../data/results/[NAME]
```

### Required Parameters

- `--analysis`: Type of analysis to perform (`kmer`, `spectrum`, `length`, `error`, `comprehensive`)
- `--output`: Output directory (automatically creates in `data/results/`)

### Optional Parameters

- Parameter ranges: `--k-range`, `--n-range`, `--length-range`, `--error-rates`
- Fixed values: `--k`, `--n`, `--length`
- Control: `--repeats`, `--step`

---

## Analysis Types

### 1. K-mer Size Analysis (`kmer`)

Tests how different k-mer sizes affect reconstruction quality.

- **Purpose**: Find optimal k-mer length for given conditions
- **Fixed**: Spectrum size (n), sequence length
- **Variable**: K-mer size (k)

### 2. Spectrum Size Analysis (`spectrum`)

Tests how the number of k-mers in the spectrum affects quality.

- **Purpose**: Determine minimum spectrum size needed
- **Fixed**: K-mer size (k), sequence length
- **Variable**: Spectrum size (n)

### 3. Sequence Length Analysis (`length`)

Tests reconstruction quality for different sequence lengths.

- **Purpose**: Understand scalability with sequence length
- **Fixed**: K-mer size (k), spectrum size (n)
- **Variable**: Sequence length

### 4. Error Tolerance Analysis (`error`)

Tests algorithm robustness against noisy data.

- **Purpose**: Determine maximum error rate the algorithm can handle
- **Fixed**: K-mer size (k), spectrum size (n), sequence length
- **Variable**: Error rate

### 5. Comprehensive Analysis (`comprehensive`)

Runs all analysis types with specified parameter ranges.

- **Purpose**: Complete parameter space exploration
- **Variable**: All parameters
- **Warning**: Can take significant time with large parameter ranges

---

## Parameter Testing Examples

### 🔬 K-mer Size Analysis

#### Basic K-mer Analysis

```bash
# Test k-mer sizes from 8 to 16 with step 2
python quality_analysis.py --analysis kmer --k-range 8 16 --n 100 --length 200 --output ../data/results/kmer_basic
```

#### Detailed K-mer Analysis

```bash
# High-resolution k-mer analysis with more repeats
python quality_analysis.py --analysis kmer \
    --k-range 6 20 --step 1 \
    --n 150 --length 250 \
    --repeats 10 \
    --output ../data/results/kmer_detailed
```

#### Quick K-mer Test

```bash
# Fast test for debugging
python quality_analysis.py --analysis kmer \
    --k-range 10 14 \
    --n 80 --length 150 \
    --repeats 3 \
    --output ../data/results/kmer_quick
```

### 📊 Spectrum Size Analysis

#### Basic Spectrum Analysis

```bash
# Test spectrum sizes from 50 to 200
python quality_analysis.py --analysis spectrum \
    --k 12 --n-range 50 200 \
    --length 200 \
    --output ../data/results/spectrum_basic
```

#### Large Spectrum Test

```bash
# Test large spectrum sizes
python quality_analysis.py --analysis spectrum \
    --k 10 --n-range 100 500 \
    --length 300 --repeats 5 \
    --output ../data/results/spectrum_large
```

#### Small Sequence Spectrum Test

```bash
# Optimize for small sequences
python quality_analysis.py --analysis spectrum \
    --k 8 --n-range 30 100 \
    --length 100 --repeats 3 \
    --output ../data/results/spectrum_small
```

### 📏 Sequence Length Analysis

#### Basic Length Analysis

```bash
# Test sequence lengths from 100-300 bp
python quality_analysis.py --analysis length \
    --k 12 --n 100 \
    --length-range 100 300 \
    --output ../data/results/length_basic
```

#### Long Sequence Analysis

```bash
# Test reconstruction of long sequences
python quality_analysis.py --analysis length \
    --k 14 --n 200 \
    --length-range 200 600 \
    --repeats 3 \
    --output ../data/results/length_long
```

#### Short Sequence Analysis

```bash
# Optimize for short sequences
python quality_analysis.py --analysis length \
    --k 8 --n 50 \
    --length-range 50 150 \
    --repeats 5 \
    --output ../data/results/length_short
```

### ⚠️ Error Tolerance Analysis

#### Basic Error Tolerance

```bash
# Test error rates from 0% to 20%
python quality_analysis.py --analysis error \
    --k 12 --n 100 --length 200 \
    --error-rates 0.0 0.05 0.1 0.15 0.2 \
    --output ../data/results/error_basic
```

#### High Error Rate Testing

```bash
# Test algorithm limits with high error rates
python quality_analysis.py --analysis error \
    --k 10 --n 150 --length 150 \
    --error-rates 0.0 0.1 0.2 0.3 0.4 0.5 \
    --repeats 10 \
    --output ../data/results/error_high
```

#### Low Error Sensitivity

```bash
# Fine-grained analysis of low error rates
python quality_analysis.py --analysis error \
    --k 14 --n 200 --length 250 \
    --error-rates 0.0 0.01 0.02 0.05 \
    --repeats 5 \
    --output ../data/results/error_low
```

### 🔬 Comprehensive Analysis

#### Full Comprehensive Analysis

```bash
# Complete parameter space exploration (takes time!)
python quality_analysis.py --analysis comprehensive \
    --k-range 8 16 --n-range 50 200 --length-range 100 300 \
    --error-rates 0.0 0.05 0.1 0.2 \
    --output ../data/results/comprehensive_full
```

#### Quick Comprehensive Analysis

```bash
# Faster comprehensive analysis
python quality_analysis.py --analysis comprehensive \
    --k-range 10 14 --n-range 80 120 --length-range 150 250 \
    --repeats 3 \
    --output ../data/results/comprehensive_quick
```

#### Focused Comprehensive Analysis

```bash
# Targeted parameter exploration
python quality_analysis.py --analysis comprehensive \
    --k-range 12 16 --step 2 \
    --n-range 100 200 --length-range 200 300 \
    --error-rates 0.0 0.1 \
    --repeats 2 \
    --output ../data/results/comprehensive_focused
```

---

## Specialized Testing Scenarios

### 🎯 Optimal Parameter Finding

#### Find Best K-mer Size

```bash
# Detailed k-mer optimization
python quality_analysis.py --analysis kmer \
    --k-range 6 20 --step 1 \
    --n 100 --length 200 \
    --repeats 5 \
    --output ../data/results/optimal_k
```

#### K-mer vs Spectrum Size Relationship

```bash
# Study interaction between k and n
python quality_analysis.py --analysis comprehensive \
    --k-range 8 16 --n-range 50 300 \
    --length 200 --repeats 3 \
    --output ../data/results/k_vs_spectrum
```

### 🔍 Edge Case Testing

#### Very Small Sequences

```bash
# Test algorithm limits with small sequences
python quality_analysis.py --analysis comprehensive \
    --k-range 4 8 --n-range 20 60 --length-range 30 80 \
    --repeats 5 \
    --output ../data/results/very_small
```

#### Large Parameter Space

```bash
# Extensive parameter exploration
python quality_analysis.py --analysis comprehensive \
    --k-range 10 20 --n-range 100 500 --length-range 200 800 \
    --repeats 2 \
    --output ../data/results/large_space
```

### ⚡ Performance Testing

#### Algorithm Limits

```bash
# Test performance with large spectra
python quality_analysis.py --analysis spectrum \
    --k 16 --n-range 300 1000 \
    --length 500 --repeats 2 \
    --output ../data/results/performance_test
```

#### Large K-mer Testing

```bash
# Test with large k-mer sizes
python quality_analysis.py --analysis kmer \
    --k-range 16 24 --step 2 \
    --n 200 --length 400 \
    --repeats 3 \
    --output ../data/results/large_kmers
```

### 🧪 Debug and Development

#### Minimal Test

```bash
# Quick verification test
python quality_analysis.py --analysis error \
    --k 10 --n 50 --length 100 \
    --error-rates 0.0 0.1 \
    --repeats 1 \
    --output ../data/results/minimal_test
```

#### Small Parameter Sweep

```bash
# Small-scale parameter testing
python quality_analysis.py --analysis kmer \
    --k-range 8 12 --step 2 \
    --n 60 --length 120 \
    --repeats 2 \
    --output ../data/results/small_sweep
```

---

## Parameter Reference

### Command Line Arguments

| Parameter        | Type     | Description                | Example                                                | Default                       |
| ---------------- | -------- | -------------------------- | ------------------------------------------------------ | ----------------------------- |
| `--analysis`     | Required | Analysis type              | `kmer`, `spectrum`, `length`, `error`, `comprehensive` | -                             |
| `--k-range`      | Optional | K-mer size range           | `8 16`                                                 | `[10, 16]`                    |
| `--n-range`      | Optional | Spectrum size range        | `50 200`                                               | `[50, 200]`                   |
| `--length-range` | Optional | Sequence length range      | `100 300`                                              | `[100, 300]`                  |
| `--error-rates`  | Optional | Error rates to test        | `0.0 0.05 0.1 0.2`                                     | `[0.0, 0.05, 0.1, 0.15, 0.2]` |
| `--k`            | Optional | Fixed k-mer size           | `12`                                                   | `12`                          |
| `--n`            | Optional | Fixed spectrum size        | `100`                                                  | `100`                         |
| `--length`       | Optional | Fixed sequence length      | `200`                                                  | `200`                         |
| `--repeats`      | Optional | Number of repeats per test | `5`                                                    | `5`                           |
| `--step`         | Optional | Step size for ranges       | `2`                                                    | `2`                           |
| `--output`       | Optional | Output directory name      | `../data/results/my_test`                              | Auto-generated timestamp      |

### Parameter Guidelines

#### K-mer Size (k)

- **Small k (4-8)**: Fast but less specific, good for short sequences
- **Medium k (10-16)**: Balanced performance, good general choice
- **Large k (18+)**: High specificity but requires longer sequences

#### Spectrum Size (n)

- **Minimum**: `sequence_length - k + 1` (for perfect coverage)
- **Recommended**: `1.5 × (sequence_length - k + 1)` (for redundancy)
- **Large n**: Better coverage but slower reconstruction

#### Sequence Length

- **Short (50-150 bp)**: Use smaller k (6-10) and appropriate n
- **Medium (150-400 bp)**: Standard parameters work well
- **Long (400+ bp)**: May need larger k and more spectrum coverage

#### Error Rates

- **0.0**: Perfect data (baseline)
- **0.05-0.1**: Realistic experimental noise
- **0.2+**: High noise conditions

---

## Best Practices

### 🎯 Planning Your Analysis

1. **Start Small**: Begin with quick tests to verify setup
2. **Single Parameter**: Test one parameter at a time first
3. **Build Up**: Gradually increase complexity
4. **Document**: Use descriptive output names

### 📊 Parameter Selection

```bash
# Step 1: Quick verification
python quality_analysis.py --analysis error --k 10 --n 50 --length 100 --error-rates 0.0 --repeats 1 --output ../data/results/verify

# Step 2: Find optimal k
python quality_analysis.py --analysis kmer --k-range 8 16 --n 100 --length 200 --repeats 3 --output ../data/results/find_k

# Step 3: Optimize spectrum size
python quality_analysis.py --analysis spectrum --k [BEST_K] --n-range 50 200 --length 200 --repeats 3 --output ../data/results/optimize_n

# Step 4: Test error tolerance
python quality_analysis.py --analysis error --k [BEST_K] --n [BEST_N] --length 200 --error-rates 0.0 0.05 0.1 0.2 --repeats 5 --output ../data/results/test_errors
```

### ⚡ Performance Optimization

#### Fast Testing

```bash
# Use fewer repeats for initial exploration
--repeats 2

# Use larger steps for range testing
--step 4

# Test smaller parameter ranges first
--k-range 10 14
```

#### Comprehensive Testing

```bash
# Use more repeats for final analysis
--repeats 10

# Use step 1 for detailed analysis
--step 1

# Test full parameter ranges
--k-range 6 20
```

### 📁 Output Organization

```bash
# Use descriptive names
--output ../data/results/kmer_optimization_20250609
--output ../data/results/error_tolerance_high_noise
--output ../data/results/spectrum_small_sequences

# Organize by date and purpose
--output ../data/results/2025_06_09/baseline_performance
--output ../data/results/2025_06_09/parameter_optimization
```

---

## Reading Results

### 📊 Result Files

Each analysis generates:

- `quality_analysis_report.txt` - Human-readable summary
- `detailed_results.csv` - Raw data for further analysis
- `analysis.log` - Execution log with timestamps

### 📈 Quick Result Check

```bash
# View summary report
cat ../data/results/your_analysis/quality_analysis_report.txt

# Check CSV data structure
head -5 ../data/results/your_analysis/detailed_results.csv

# Check for errors in log
grep -i error ../data/results/your_analysis/analysis.log
```

### 📋 Understanding Output

#### Report Structure

```
=== DNA Sequence Reconstruction Quality Analysis Report ===
Analysis Date: 2025-06-09 16:54:21
Total Tests Performed: X

=== Summary Statistics ===
Overall Performance:
  Average Accuracy: XX.XX% ± X.XX%
  Average Coverage: XX.XX% ± X.XX%
  Average Runtime: X.XXXXs ± X.XXXXs
  Success Rate: XXX.XX%

Best Performance Configurations:
  Highest Accuracy: XX.XX% (k=XX, n=XXX, length=XXX)
  Fastest Runtime: X.XXXXs (k=XX, n=XXX, length=XXX)

=== Recommendations for High Quality ===
[Specific recommendations based on results]
```

#### CSV Columns

```
k,n,seq_length,error_rate,original_length,reconstructed_length,coverage,accuracy,edit_distance,runtime,is_valid,success,repeat
```

---

## Troubleshooting

### ❌ Common Issues

#### Command Not Found

```bash
# Make sure you're in the right directory
cd scripts
ls -la  # Should see quality_analysis.py
```

#### No Results Generated

```bash
# Check if analysis completed
cat ../data/results/your_analysis/analysis.log | tail -10

# Check for permission issues
ls -la ../data/results/
```

#### Infinite Loop (Fixed!)

The backtracking infinite loop issue has been resolved. The algorithm now:

- ✅ Tracks failed paths to avoid repeating them
- ✅ Has timeout protection (30 seconds max)
- ✅ Limits backtrack attempts (max 10)
- ✅ Includes iteration counters

#### Memory Issues

```bash
# Reduce parameter ranges
--repeats 2
--step 4

# Test smaller sequences first
--length-range 50 150
```

### 🔧 Debug Commands

#### Minimal Test

```bash
python quality_analysis.py --analysis error --k 8 --n 30 --length 50 --error-rates 0.0 --repeats 1 --output ../data/results/debug
```

#### Verbose Logging

```bash
# Check the analysis.log for detailed information
tail -f ../data/results/your_analysis/analysis.log
```

### 📞 Getting Help

```bash
# View all available options
python quality_analysis.py --help

# View examples
python quality_analysis.py --help | grep -A 20 "Examples:"
```

---

## Advanced Usage

### 🔄 Batch Processing

Create a script for multiple analyses:

```bash
#!/bin/bash
# batch_analysis.sh

# Set base output directory
BASE_DIR="../data/results/batch_$(date +%Y%m%d_%H%M%S)"

# Run multiple analyses
python quality_analysis.py --analysis kmer --k-range 8 16 --n 100 --length 200 --output "${BASE_DIR}/kmer_analysis"
python quality_analysis.py --analysis spectrum --k 12 --n-range 50 200 --length 200 --output "${BASE_DIR}/spectrum_analysis"
python quality_analysis.py --analysis error --k 12 --n 100 --length 200 --error-rates 0.0 0.05 0.1 0.2 --output "${BASE_DIR}/error_analysis"

echo "Batch analysis completed. Results in: $BASE_DIR"
```

### 📊 Data Analysis

After generating CSV files, you can analyze them with:

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load results
df = pd.read_csv('../data/results/your_analysis/detailed_results.csv')

# Plot accuracy vs k-mer size
df.groupby('k')['accuracy'].mean().plot(kind='bar')
plt.title('Accuracy vs K-mer Size')
plt.ylabel('Accuracy (%)')
plt.show()

# Find best parameters
best = df.loc[df['accuracy'].idxmax()]
print(f"Best configuration: k={best['k']}, n={best['n']}, accuracy={best['accuracy']:.2f}%")
```

---

## Summary

This guide provides comprehensive examples for testing all aspects of the DNA sequence reconstruction algorithm. Start with quick tests, then progressively explore the parameter space to find optimal configurations for your specific use case.

**Key Takeaways:**

- Always start with small, quick tests
- Use descriptive output names for organization
- Test one parameter at a time before comprehensive analysis
- Check results in both report and CSV format
- The infinite loop issue has been fixed - all tests should complete
- All results are automatically saved to `data/results/` for easy access

Happy parameter testing! 🧬🔬
