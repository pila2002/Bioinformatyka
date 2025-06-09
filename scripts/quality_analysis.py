#!/usr/bin/env python3
"""
Quality Analysis Script for DNA Sequence Reconstruction

This script performs comprehensive analysis of reconstruction quality
against various parameters to help understand how to achieve high-quality results.
"""

import argparse
import sys
import os
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from typing import List, Dict, Tuple
import json

# Add src directory to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from algorithms.classic_sbh import ClassicSBH
from utils.sequence_generator import SequenceGenerator
from utils.spectrum_generator import SpectrumGenerator
from utils.sequence_validator import SequenceValidator
from utils.logger import Logger


class QualityAnalyzer:
    """Comprehensive quality analyzer for DNA sequence reconstruction."""
    
    def __init__(self, output_dir: str, logger: Logger):
        self.output_dir = output_dir
        self.logger = logger
        self.results = []
        
    def run_single_test(self, k: int, n: int, seq_length: int, error_rate: float = 0.0) -> Dict:
        """Run a single reconstruction test and return quality metrics."""
        
        # Generate test sequence
        generator = SequenceGenerator()
        original_sequence = generator.generate_random_sequence(seq_length)
        
        # Generate spectrum
        spectrum_gen = SpectrumGenerator()
        full_spectrum = spectrum_gen.generate_spectrum(original_sequence, k)
        
        # Sample n k-mers from the full spectrum if needed
        if n < len(full_spectrum):
            import random
            spectrum = random.sample(list(full_spectrum), n)
        else:
            spectrum = list(full_spectrum)
        
        # Add errors if specified
        if error_rate > 0:
            spectrum = self._add_spectrum_errors(spectrum, error_rate)
        
        # Perform reconstruction
        algorithm = ClassicSBH()
        validator = SequenceValidator()
        
        start_time = time.time()
        try:
            reconstructed_sequence = algorithm.reconstruct(spectrum, seq_length, k)
            runtime = time.time() - start_time
            
            # Calculate quality metrics
            is_valid, coverage = validator.validate_reconstruction(reconstructed_sequence, spectrum, k)
            accuracy = self._calculate_accuracy(original_sequence, reconstructed_sequence)
            edit_distance = self._calculate_edit_distance(original_sequence, reconstructed_sequence)
            
            return {
                'k': k,
                'n': n,
                'seq_length': seq_length,
                'error_rate': error_rate,
                'original_length': len(original_sequence),
                'reconstructed_length': len(reconstructed_sequence) if reconstructed_sequence else 0,
                'coverage': coverage,
                'accuracy': accuracy,
                'edit_distance': edit_distance,
                'runtime': runtime,
                'is_valid': is_valid,
                'success': reconstructed_sequence is not None
            }
            
        except Exception as e:
            runtime = time.time() - start_time
            self.logger.log(f"Reconstruction failed: {str(e)}")
            return {
                'k': k,
                'n': n,
                'seq_length': seq_length,
                'error_rate': error_rate,
                'original_length': seq_length,
                'reconstructed_length': 0,
                'coverage': 0.0,
                'accuracy': 0.0,
                'edit_distance': seq_length,
                'runtime': runtime,
                'is_valid': False,
                'success': False
            }
    
    def _add_spectrum_errors(self, spectrum: List[str], error_rate: float) -> List[str]:
        """Add random errors to spectrum."""
        import random
        
        if error_rate == 0:
            return spectrum
            
        corrupted_spectrum = spectrum.copy()
        num_errors = int(len(spectrum) * error_rate)
        
        for _ in range(num_errors):
            if corrupted_spectrum:
                # Randomly remove or modify k-mers
                if random.random() < 0.5:  # Remove k-mer
                    idx = random.randint(0, len(corrupted_spectrum) - 1)
                    corrupted_spectrum.pop(idx)
                else:  # Modify k-mer
                    idx = random.randint(0, len(corrupted_spectrum) - 1)
                    kmer = corrupted_spectrum[idx]
                    pos = random.randint(0, len(kmer) - 1)
                    nucleotides = ['A', 'T', 'G', 'C']
                    new_nucleotide = random.choice([n for n in nucleotides if n != kmer[pos]])
                    corrupted_spectrum[idx] = kmer[:pos] + new_nucleotide + kmer[pos+1:]
        
        return corrupted_spectrum
    
    def _calculate_accuracy(self, original: str, reconstructed: str) -> float:
        """Calculate sequence accuracy (percentage of correct positions)."""
        if not reconstructed or not original:
            return 0.0
            
        min_len = min(len(original), len(reconstructed))
        if min_len == 0:
            return 0.0
            
        matches = sum(1 for i in range(min_len) if original[i] == reconstructed[i])
        return (matches / min_len) * 100.0
    
    def _calculate_edit_distance(self, s1: str, s2: str) -> int:
        """Calculate Levenshtein edit distance between two strings."""
        if not s1:
            return len(s2) if s2 else 0
        if not s2:
            return len(s1)
            
        # Create matrix
        matrix = [[0] * (len(s2) + 1) for _ in range(len(s1) + 1)]
        
        # Initialize first row and column
        for i in range(len(s1) + 1):
            matrix[i][0] = i
        for j in range(len(s2) + 1):
            matrix[0][j] = j
        
        # Fill matrix
        for i in range(1, len(s1) + 1):
            for j in range(1, len(s2) + 1):
                if s1[i-1] == s2[j-1]:
                    matrix[i][j] = matrix[i-1][j-1]
                else:
                    matrix[i][j] = min(
                        matrix[i-1][j] + 1,      # deletion
                        matrix[i][j-1] + 1,      # insertion
                        matrix[i-1][j-1] + 1     # substitution
                    )
        
        return matrix[len(s1)][len(s2)]
    
    def analyze_k_mer_size_impact(self, k_range: range, n: int, seq_length: int, 
                                  repeats: int = 5) -> pd.DataFrame:
        """Analyze impact of k-mer size on reconstruction quality."""
        self.logger.log("=== Analyzing K-mer Size Impact ===")
        
        results = []
        for k in k_range:
            self.logger.log(f"Testing k-mer size: {k}")
            for repeat in range(repeats):
                result = self.run_single_test(k, n, seq_length)
                result['repeat'] = repeat
                results.append(result)
                
        df = pd.DataFrame(results)
        self.results.extend(results)
        return df
    
    def analyze_spectrum_size_impact(self, k: int, n_range: range, seq_length: int, 
                                    repeats: int = 5) -> pd.DataFrame:
        """Analyze impact of spectrum size on reconstruction quality."""
        self.logger.log("=== Analyzing Spectrum Size Impact ===")
        
        results = []
        for n in n_range:
            self.logger.log(f"Testing spectrum size: {n}")
            for repeat in range(repeats):
                result = self.run_single_test(k, n, seq_length)
                result['repeat'] = repeat
                results.append(result)
                
        df = pd.DataFrame(results)
        self.results.extend(results)
        return df
    
    def analyze_sequence_length_impact(self, k: int, n: int, length_range: range, 
                                      repeats: int = 5) -> pd.DataFrame:
        """Analyze impact of sequence length on reconstruction quality."""
        self.logger.log("=== Analyzing Sequence Length Impact ===")
        
        results = []
        for seq_len in length_range:
            self.logger.log(f"Testing sequence length: {seq_len}")
            for repeat in range(repeats):
                result = self.run_single_test(k, n, seq_len)
                result['repeat'] = repeat
                results.append(result)
                
        df = pd.DataFrame(results)
        self.results.extend(results)
        return df
    
    def analyze_error_tolerance(self, k: int, n: int, seq_length: int, 
                               error_rates: List[float], repeats: int = 5) -> pd.DataFrame:
        """Analyze reconstruction quality under different error rates."""
        self.logger.log("=== Analyzing Error Tolerance ===")
        
        results = []
        for error_rate in error_rates:
            self.logger.log(f"Testing error rate: {error_rate:.2%}")
            for repeat in range(repeats):
                result = self.run_single_test(k, n, seq_length, error_rate)
                result['repeat'] = repeat
                results.append(result)
                
        df = pd.DataFrame(results)
        self.results.extend(results)
        return df
    
    def generate_comprehensive_report(self, results_df: pd.DataFrame):
        """Generate comprehensive analysis report."""
        self.logger.log("=== Generating Comprehensive Report ===")
        
        # Create visualizations
        self._create_quality_plots(results_df)
        
        # Generate summary statistics
        summary = self._generate_summary_stats(results_df)
        
        # Save detailed results
        results_df.to_csv(os.path.join(self.output_dir, 'detailed_results.csv'), index=False)
        
        # Save summary report
        with open(os.path.join(self.output_dir, 'quality_analysis_report.txt'), 'w') as f:
            f.write("=== DNA Sequence Reconstruction Quality Analysis Report ===\n\n")
            f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total Tests Performed: {len(results_df)}\n\n")
            
            f.write("=== Summary Statistics ===\n")
            f.write(summary)
            
            f.write("\n=== Recommendations for High Quality ===\n")
            recommendations = self._generate_recommendations(results_df)
            f.write(recommendations)
        
        self.logger.log(f"Report saved to: {self.output_dir}")
    
    def _create_quality_plots(self, df: pd.DataFrame):
        """Create comprehensive quality visualization plots."""
        try:
            plt.style.use('seaborn-v0_8')
        except:
            # Fallback to default style if seaborn style is not available
            plt.style.use('default')
        
        # Quality vs K-mer size
        if 'k' in df.columns and df['k'].nunique() > 1:
            try:
                plt.figure(figsize=(12, 8))
                
                plt.subplot(2, 2, 1)
                k_stats = df.groupby('k').agg({
                    'accuracy': ['mean', 'std'],
                    'coverage': ['mean', 'std'],
                    'success': 'mean'
                }).round(3)
                
                plt.plot(k_stats.index, k_stats[('accuracy', 'mean')], 'o-', label='Accuracy')
                plt.fill_between(k_stats.index, 
                               k_stats[('accuracy', 'mean')] - k_stats[('accuracy', 'std')],
                               k_stats[('accuracy', 'mean')] + k_stats[('accuracy', 'std')],
                               alpha=0.3)
                plt.xlabel('K-mer Size')
                plt.ylabel('Accuracy (%)')
                plt.title('Accuracy vs K-mer Size')
                plt.legend()
                plt.grid(True, alpha=0.3)
                
                plt.subplot(2, 2, 2)
                plt.plot(k_stats.index, k_stats[('coverage', 'mean')], 'o-', color='orange', label='Coverage')
                plt.fill_between(k_stats.index, 
                               k_stats[('coverage', 'mean')] - k_stats[('coverage', 'std')],
                               k_stats[('coverage', 'mean')] + k_stats[('coverage', 'std')],
                               alpha=0.3, color='orange')
                plt.xlabel('K-mer Size')
                plt.ylabel('Coverage (%)')
                plt.title('Coverage vs K-mer Size')
                plt.legend()
                plt.grid(True, alpha=0.3)
                
                plt.subplot(2, 2, 3)
                plt.plot(k_stats.index, k_stats[('success', 'mean')] * 100, 'o-', color='green', label='Success Rate')
                plt.xlabel('K-mer Size')
                plt.ylabel('Success Rate (%)')
                plt.title('Success Rate vs K-mer Size')
                plt.legend()
                plt.grid(True, alpha=0.3)
                
                plt.subplot(2, 2, 4)
                runtime_stats = df.groupby('k')['runtime'].agg(['mean', 'std'])
                plt.plot(runtime_stats.index, runtime_stats['mean'], 'o-', color='red', label='Runtime')
                plt.fill_between(runtime_stats.index, 
                               runtime_stats['mean'] - runtime_stats['std'],
                               runtime_stats['mean'] + runtime_stats['std'],
                               alpha=0.3, color='red')
                plt.xlabel('K-mer Size')
                plt.ylabel('Runtime (seconds)')
                plt.title('Runtime vs K-mer Size')
                plt.legend()
                plt.grid(True, alpha=0.3)
                
                plt.tight_layout()
                plt.savefig(os.path.join(self.output_dir, 'quality_vs_kmer_size.png'), dpi=300, bbox_inches='tight')
                plt.close()
            except Exception as e:
                self.logger.log(f"Error creating k-mer size plot: {str(e)}")
                plt.close()
        
        # Quality vs Spectrum size
        if 'n' in df.columns and df['n'].nunique() > 1:
            plt.figure(figsize=(12, 6))
            
            n_stats = df.groupby('n').agg({
                'accuracy': ['mean', 'std'],
                'coverage': ['mean', 'std'],
                'runtime': ['mean', 'std']
            }).round(3)
            
            plt.subplot(1, 2, 1)
            plt.plot(n_stats.index, n_stats[('accuracy', 'mean')], 'o-', label='Accuracy')
            plt.plot(n_stats.index, n_stats[('coverage', 'mean')], 'o-', label='Coverage')
            plt.xlabel('Spectrum Size (n)')
            plt.ylabel('Quality (%)')
            plt.title('Quality vs Spectrum Size')
            plt.legend()
            plt.grid(True, alpha=0.3)
            
            plt.subplot(1, 2, 2)
            plt.plot(n_stats.index, n_stats[('runtime', 'mean')], 'o-', color='red', label='Runtime')
            plt.xlabel('Spectrum Size (n)')
            plt.ylabel('Runtime (seconds)')
            plt.title('Runtime vs Spectrum Size')
            plt.legend()
            plt.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, 'quality_vs_spectrum_size.png'), dpi=300, bbox_inches='tight')
            plt.close()
        
        # Error tolerance analysis
        if 'error_rate' in df.columns and df['error_rate'].nunique() > 1:
            plt.figure(figsize=(10, 6))
            
            error_stats = df.groupby('error_rate').agg({
                'accuracy': ['mean', 'std'],
                'success': 'mean'
            }).round(3)
            
            plt.subplot(1, 2, 1)
            plt.plot(error_stats.index * 100, error_stats[('accuracy', 'mean')], 'o-', label='Accuracy')
            plt.xlabel('Error Rate (%)')
            plt.ylabel('Accuracy (%)')
            plt.title('Accuracy vs Error Rate')
            plt.legend()
            plt.grid(True, alpha=0.3)
            
            plt.subplot(1, 2, 2)
            plt.plot(error_stats.index * 100, error_stats[('success', 'mean')] * 100, 'o-', color='green', label='Success Rate')
            plt.xlabel('Error Rate (%)')
            plt.ylabel('Success Rate (%)')
            plt.title('Success Rate vs Error Rate')
            plt.legend()
            plt.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, 'error_tolerance.png'), dpi=300, bbox_inches='tight')
            plt.close()
    
    def _generate_summary_stats(self, df: pd.DataFrame) -> str:
        """Generate summary statistics report."""
        report = []
        
        if len(df) == 0:
            report.append("No data available for analysis.")
            return "\n".join(report)
        
        # Calculate statistics safely
        accuracy_mean = df['accuracy'].mean()
        accuracy_std = df['accuracy'].std()
        coverage_mean = df['coverage'].mean()
        coverage_std = df['coverage'].std()
        runtime_mean = df['runtime'].mean()
        runtime_std = df['runtime'].std()
        success_rate = df['success'].mean()
        
        report.append("Overall Performance:")
        report.append(f"  Average Accuracy: {accuracy_mean:.2f}% ± {accuracy_std:.2f}%")
        report.append(f"  Average Coverage: {coverage_mean:.2f}% ± {coverage_std:.2f}%")
        report.append(f"  Average Runtime: {runtime_mean:.4f}s ± {runtime_std:.4f}s")
        report.append(f"  Success Rate: {success_rate * 100:.2f}%")
        report.append("")
        
        # Best performance analysis
        if len(df) > 0 and df['accuracy'].max() > 0:
            best_accuracy = df.loc[df['accuracy'].idxmax()]
            best_runtime = df.loc[df['runtime'].idxmin()]
            
            report.append("Best Performance Configurations:")
            report.append(f"  Highest Accuracy: {best_accuracy['accuracy']:.2f}% (k={best_accuracy['k']}, n={best_accuracy['n']}, length={best_accuracy['seq_length']})")
            report.append(f"  Fastest Runtime: {best_runtime['runtime']:.4f}s (k={best_runtime['k']}, n={best_runtime['n']}, length={best_runtime['seq_length']})")
            report.append("")
        else:
            report.append("All reconstruction attempts failed.")
            report.append("")
        
        return "\n".join(report)
    
    def _generate_recommendations(self, df: pd.DataFrame) -> str:
        """Generate recommendations for achieving high quality."""
        recommendations = []
        
        # Analyze k-mer size impact
        if 'k' in df.columns and df['k'].nunique() > 1:
            k_quality = df.groupby('k')['accuracy'].mean()
            best_k = k_quality.idxmax()
            recommendations.append(f"1. Optimal K-mer Size: Use k={best_k} for best accuracy ({k_quality[best_k]:.2f}%)")
        
        # Analyze spectrum size impact
        if 'n' in df.columns and df['n'].nunique() > 1:
            n_quality = df.groupby('n')['accuracy'].mean()
            best_n = n_quality.idxmax()
            recommendations.append(f"2. Optimal Spectrum Size: Use n={best_n} k-mers for best results")
        
        # Error tolerance analysis
        if 'error_rate' in df.columns and df['error_rate'].nunique() > 1:
            error_tolerance = df[df['success'] == True]['error_rate'].max()
            recommendations.append(f"3. Error Tolerance: Algorithm can handle up to {error_tolerance:.2%} error rate")
        
        # General recommendations
        recommendations.append("4. General Guidelines:")
        recommendations.append("   - Higher k-mer sizes generally provide better specificity but require longer sequences")
        recommendations.append("   - Ensure spectrum size (n) is at least equal to sequence length - k + 1")
        recommendations.append("   - For noisy data, consider using error correction or multiple reconstructions")
        recommendations.append("   - Monitor coverage percentage - aim for 100% coverage for best results")
        
        return "\n".join(recommendations)


def main():
    parser = argparse.ArgumentParser(
        description="Comprehensive Quality Analysis for DNA Sequence Reconstruction",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic k-mer size analysis
  python quality_analysis.py --analysis kmer --k-range 8 16 --n 100 --length 200

  # Spectrum size analysis
  python quality_analysis.py --analysis spectrum --k 12 --n-range 50 200 --length 150

  # Error tolerance analysis
  python quality_analysis.py --analysis error --k 10 --n 100 --length 100 --error-rates 0.0 0.1 0.2

  # Comprehensive analysis (all parameters)
  python quality_analysis.py --analysis comprehensive --k-range 8 16 --n-range 50 200 --length-range 100 300
        """
    )
    
    parser.add_argument('--analysis', choices=['kmer', 'spectrum', 'length', 'error', 'comprehensive'], 
                       required=True, help='Type of analysis to perform')
    
    # Parameter ranges
    parser.add_argument('--k-range', nargs=2, type=int, default=[10, 16], 
                       help='K-mer size range (start, end)')
    parser.add_argument('--n-range', nargs=2, type=int, default=[50, 200], 
                       help='Spectrum size range (start, end)')
    parser.add_argument('--length-range', nargs=2, type=int, default=[100, 300], 
                       help='Sequence length range (start, end)')
    parser.add_argument('--error-rates', nargs='+', type=float, 
                       default=[0.0, 0.05, 0.1, 0.15, 0.2], 
                       help='Error rates to test')
    
    # Fixed parameters for single-parameter analysis
    parser.add_argument('--k', type=int, default=12, help='Fixed k-mer size')
    parser.add_argument('--n', type=int, default=100, help='Fixed spectrum size')
    parser.add_argument('--length', type=int, default=200, help='Fixed sequence length')
    
    # Analysis parameters
    parser.add_argument('--repeats', type=int, default=5, 
                       help='Number of repeats per configuration')
    parser.add_argument('--step', type=int, default=2, 
                       help='Step size for parameter ranges')
    
    # Output parameters
    parser.add_argument('--output', type=str, default=None, 
                       help='Output directory (default: auto-generated)')
    
    args = parser.parse_args()
    
    # Create output directory
    if args.output is None:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_dir = f'quality_analysis_{timestamp}'
    else:
        output_dir = args.output
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize logger and analyzer
    logger = Logger(os.path.join(output_dir, 'analysis.log'))
    analyzer = QualityAnalyzer(output_dir, logger)
    
    logger.log(f"Starting quality analysis: {args.analysis}")
    logger.log(f"Output directory: {output_dir}")
    
    all_results = []
    
    try:
        if args.analysis == 'kmer':
            k_range = range(args.k_range[0], args.k_range[1] + 1, args.step)
            df = analyzer.analyze_k_mer_size_impact(k_range, args.n, args.length, args.repeats)
            all_results.append(df)
            
        elif args.analysis == 'spectrum':
            n_range = range(args.n_range[0], args.n_range[1] + 1, 
                           max(1, (args.n_range[1] - args.n_range[0]) // 10))
            df = analyzer.analyze_spectrum_size_impact(args.k, n_range, args.length, args.repeats)
            all_results.append(df)
            
        elif args.analysis == 'length':
            length_range = range(args.length_range[0], args.length_range[1] + 1, 
                               max(1, (args.length_range[1] - args.length_range[0]) // 10))
            df = analyzer.analyze_sequence_length_impact(args.k, args.n, length_range, args.repeats)
            all_results.append(df)
            
        elif args.analysis == 'error':
            df = analyzer.analyze_error_tolerance(args.k, args.n, args.length, 
                                                 args.error_rates, args.repeats)
            all_results.append(df)
            
        elif args.analysis == 'comprehensive':
            # Run all analyses
            k_range = range(args.k_range[0], args.k_range[1] + 1, args.step)
            n_range = range(args.n_range[0], args.n_range[1] + 1, 
                           max(1, (args.n_range[1] - args.n_range[0]) // 5))
            length_range = range(args.length_range[0], args.length_range[1] + 1, 
                               max(1, (args.length_range[1] - args.length_range[0]) // 5))
            
            df1 = analyzer.analyze_k_mer_size_impact(k_range, args.n, args.length, args.repeats)
            df2 = analyzer.analyze_spectrum_size_impact(args.k, n_range, args.length, args.repeats)
            df3 = analyzer.analyze_sequence_length_impact(args.k, args.n, length_range, args.repeats)
            df4 = analyzer.analyze_error_tolerance(args.k, args.n, args.length, 
                                                  args.error_rates, args.repeats)
            
            all_results.extend([df1, df2, df3, df4])
        
        # Combine all results
        if all_results:
            combined_df = pd.concat(all_results, ignore_index=True)
            analyzer.generate_comprehensive_report(combined_df)
            
            logger.log(f"Analysis completed successfully!")
            logger.log(f"Results saved to: {output_dir}")
            print(f"Quality analysis completed! Results saved to: {output_dir}")
        
    except Exception as e:
        logger.log(f"Error during analysis: {str(e)}")
        print(f"Error: {str(e)}")
        return 1
    
    return 0


if __name__ == '__main__':
    sys.exit(main()) 