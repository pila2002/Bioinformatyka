# Documentation Index

Welcome to the DNA Sequence Reconstruction project documentation. This directory contains comprehensive guides and documentation for using and understanding the SBH (Sequencing by Hybridization) algorithms.

## 📚 Available Documentation

### 🎯 [Parameter Testing Guide](parameter_testing_guide.md)

**Complete guide for testing different algorithm parameters**

- Command-line usage examples
- Analysis types (k-mer, spectrum, length, error, comprehensive)
- Parameter optimization strategies
- Best practices and troubleshooting
- Advanced usage and batch processing

**Quick Start:**

```bash
cd scripts
python quality_analysis.py --analysis error --k 10 --n 50 --length 100 --error-rates 0.0 0.1 --repeats 1 --output ../data/results/quick_test
```

### 📊 [Report Directory](report/)

Contains generated analysis reports and documentation from previous runs.

## 🚀 Getting Started

1. **First Time Setup**: Check the main project README for environment setup
2. **Quick Test**: Use the parameter testing guide to run your first analysis
3. **Parameter Optimization**: Follow the step-by-step optimization strategy
4. **Results Analysis**: Learn how to interpret and analyze your results

## 🔗 Related Files

- **Main Project**: [`../project_outline.md`](../project_outline.md) - Overall project structure
- **Data Organization**: [`../data/README.md`](../data/README.md) - Data folder organization
- **Source Code**: [`../src/`](../src/) - Algorithm implementations
- **Scripts**: [`../scripts/`](../scripts/) - Analysis and testing scripts

## 📈 Common Use Cases

### Finding Optimal Parameters

1. Start with quick k-mer analysis: `--analysis kmer --k-range 8 16`
2. Optimize spectrum size: `--analysis spectrum --n-range 50 200`
3. Test error tolerance: `--analysis error --error-rates 0.0 0.05 0.1 0.2`

### Performance Testing

1. Algorithm limits: Test with large parameter ranges
2. Speed optimization: Use fewer repeats for initial exploration
3. Memory usage: Start with smaller sequences and ranges

### Research and Development

1. Comprehensive analysis: `--analysis comprehensive` for full parameter space
2. Batch processing: Use the provided batch script templates
3. Data analysis: Import CSV results into pandas/R for further analysis

## 🛠️ Recent Updates

### 2025-06-09

- ✅ **Fixed infinite loop issue** in backtracking algorithm
- ✅ **Added comprehensive parameter testing guide** with 50+ examples
- ✅ **Organized data folder structure** for better result management
- ✅ **Enhanced error handling** and timeout protection
- ✅ **Improved logging** for better debugging

## 🆘 Need Help?

1. **Parameter Testing**: See [parameter_testing_guide.md](parameter_testing_guide.md)
2. **Algorithm Issues**: Check the troubleshooting section in the parameter guide
3. **Data Analysis**: Use the CSV analysis examples in the parameter guide
4. **Performance**: Try the quick test examples first

## 📝 Contributing

When adding new documentation:

1. Follow the existing markdown structure
2. Include practical examples with code blocks
3. Add cross-references to related documentation
4. Update this index file with new content

---

**Happy researching!** 🧬🔬
