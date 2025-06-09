import argparse
import time
from algorithms.classic_sbh import ClassicSBH
from utils.sequence_generator import SequenceGenerator
from utils.spectrum_generator import SpectrumGenerator
from utils.sequence_validator import SequenceValidator
from utils.logger import Logger
from utils.visualization import BenchmarkVisualizer
from typing import List, Dict
import json
import os
import pandas as pd
from datetime import datetime

def run_benchmark(k: int, n: int, num_tests: int, logger: Logger) -> List[Dict]:
    """Run benchmark tests and return results."""
    results = []
    
    # Initialize components
    sequence_generator = SequenceGenerator()
    spectrum_generator = SpectrumGenerator()
    sequence_validator = SequenceValidator()
    sbh = ClassicSBH()
    
    for i in range(num_tests):
        logger.log(f"\nTest case {i+1}/{num_tests}")
        logger.log(f"Parameters: k={k}, n={n}")
        
        # Generate random sequence
        original_sequence = sequence_generator.generate_random_sequence(n + k - 1)
        logger.log(f"Original sequence: {original_sequence}")

        # Generate spectrum
        spectrum = spectrum_generator.generate_spectrum(original_sequence, k)
        logger.log(f"Generated spectrum with {len(spectrum)} k-mers")

        # Measure runtime
        start_time = time.time()
        reconstructed_sequence = sbh.reconstruct(list(spectrum), len(original_sequence), k)
        runtime = time.time() - start_time
        
        logger.log(f"Reconstructed sequence: {reconstructed_sequence}")

        # Validate reconstruction
        is_valid, coverage = sequence_validator.validate_reconstruction(
            reconstructed_sequence, spectrum, k
        )
        logger.log(f"Valid: {is_valid}")
        logger.log(f"Spectrum coverage in reconstruction: {coverage:.2f}%")
        
        # Store results
        results.append({
            'k': k,
            'n': n,
            'sequence_length': len(original_sequence),
            'coverage': coverage,
            'runtime': runtime,
            'is_valid': is_valid
        })
    
    return results

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='DNA Sequence Reconstruction')
    parser.add_argument('--k', type=int, default=10, help='k-mer size')
    parser.add_argument('--n', type=int, default=100, help='number of k-mers')
    parser.add_argument('--tests', type=int, default=5, help='number of test cases')
    parser.add_argument('--benchmark', action='store_true', help='run benchmark tests')
    args = parser.parse_args()
    
    # Create results directory
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    results_dir = os.path.join('benchmark_results', timestamp)
    os.makedirs(results_dir, exist_ok=True)
    
    # Initialize logger
    logger = Logger(os.path.join(results_dir, 'benchmark.log'))
    
    if args.benchmark:
        # Run benchmark with multiple k values
        k_values = [8, 10, 12, 14, 16]
        all_results = []
        
        for k in k_values:
            logger.log(f"\n=== Running benchmark for k={k} ===")
            results = run_benchmark(k, args.n, args.tests, logger)
            all_results.extend(results)
        
        # Save raw results in both JSON and CSV formats
        with open(os.path.join(results_dir, 'raw_results.json'), 'w') as f:
            json.dump(all_results, f, indent=2)
        
        # Convert to DataFrame and save as CSV
        df = pd.DataFrame(all_results)
        df.to_csv(os.path.join(results_dir, 'raw_results.csv'), index=False)
        
        # Generate visualizations
        visualizer = BenchmarkVisualizer(os.path.join(results_dir, 'plots'))
        visualizer.plot_all_metrics(all_results)
        
        logger.log(f"\nBenchmark results saved to {results_dir}")
    else:
        # Run single test
        results = run_benchmark(args.k, args.n, args.tests, logger)
        
        # Save results in both JSON and CSV formats
        with open(os.path.join(results_dir, 'results.json'), 'w') as f:
            json.dump(results, f, indent=2)
        
        # Convert to DataFrame and save as CSV
        df = pd.DataFrame(results)
        df.to_csv(os.path.join(results_dir, 'results.csv'), index=False)
        
        # Generate visualizations
        visualizer = BenchmarkVisualizer(os.path.join(results_dir, 'plots'))
        visualizer.plot_all_metrics(results)
        
        logger.log(f"\nResults saved to {results_dir}")

if __name__ == '__main__':
    main() 