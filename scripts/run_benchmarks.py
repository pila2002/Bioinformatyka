import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from pathlib import Path
import logging

# Add src to Python path
sys.path.append(str(Path(__file__).parent.parent))

from src.benchmarking.benchmark_sbh import SBHBenchmark, BenchmarkParameters

def setup_logging(output_dir: Path) -> None:
    """Set up logging configuration"""
    log_file = output_dir / 'benchmark.log'
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

def run_length_vs_accuracy_benchmark():
    """Benchmark sequence length vs reconstruction accuracy"""
    benchmark = SBHBenchmark()
    results = []
    
    # Test different sequence lengths
    sequence_lengths = [100, 200, 300, 400, 500]
    k_mer_lengths = [6, 8, 10]
    
    for seq_len in sequence_lengths:
        for k in k_mer_lengths:
            params = BenchmarkParameters(
                sequence_length=seq_len,
                k_mer_length=k,
                negative_error_rate=0.0,
                positive_error_rate=0.0,
                num_trials=10
            )
            
            logging.info(f"Running length vs accuracy benchmark: length={seq_len}, k={k}")
            result = benchmark.run_benchmark(params)
            logging.info(f"Results: accuracy={result.avg_reconstruction_accuracy:.3f} ± {result.std_reconstruction_accuracy:.3f}, "
                        f"time={result.avg_execution_time:.3f}s ± {result.std_execution_time:.3f}s, "
                        f"success_rate={result.successful_reconstructions/params.num_trials:.2%}")
            
            results.append({
                'sequence_length': seq_len,
                'k_mer_length': k,
                'accuracy': result.avg_reconstruction_accuracy,
                'time': result.avg_execution_time,
                'accuracy_std': result.std_reconstruction_accuracy,
                'time_std': result.std_execution_time,
                'success_rate': result.successful_reconstructions / params.num_trials
            })
    
    return pd.DataFrame(results)

def run_error_rate_benchmark():
    """Benchmark error rates vs reconstruction accuracy"""
    benchmark = SBHBenchmark()
    results = []
    
    # Test different error rates
    error_rates = [0.0, 0.05, 0.1, 0.15, 0.2]
    sequence_length = 300
    k_mer_length = 8
    
    for neg_error in error_rates:
        for pos_error in error_rates:
            params = BenchmarkParameters(
                sequence_length=sequence_length,
                k_mer_length=k_mer_length,
                negative_error_rate=neg_error,
                positive_error_rate=pos_error,
                num_trials=10
            )
            
            logging.info(f"Running error rate benchmark: neg_error={neg_error}, pos_error={pos_error}")
            result = benchmark.run_benchmark(params)
            logging.info(f"Results: accuracy={result.avg_reconstruction_accuracy:.3f} ± {result.std_reconstruction_accuracy:.3f}, "
                        f"time={result.avg_execution_time:.3f}s ± {result.std_execution_time:.3f}s, "
                        f"success_rate={result.successful_reconstructions/params.num_trials:.2%}")
            
            results.append({
                'negative_error_rate': neg_error,
                'positive_error_rate': pos_error,
                'accuracy': result.avg_reconstruction_accuracy,
                'time': result.avg_execution_time,
                'accuracy_std': result.std_reconstruction_accuracy,
                'time_std': result.std_execution_time,
                'success_rate': result.successful_reconstructions / params.num_trials
            })
    
    return pd.DataFrame(results)

def plot_length_vs_accuracy(df: pd.DataFrame, output_dir: Path):
    """Plot sequence length vs accuracy results"""
    plt.figure(figsize=(12, 6))
    
    for k in df['k_mer_length'].unique():
        data = df[df['k_mer_length'] == k]
        plt.errorbar(
            data['sequence_length'],
            data['accuracy'],
            yerr=data['accuracy_std'],
            label=f'k={k}',
            marker='o'
        )
    
    plt.xlabel('Sequence Length')
    plt.ylabel('Reconstruction Accuracy')
    plt.title('Sequence Length vs Reconstruction Accuracy')
    plt.legend()
    plt.grid(True)
    plt.savefig(output_dir / 'length_vs_accuracy.png')
    plt.close()
    
    logging.info(f"Length vs accuracy plot saved to {output_dir / 'length_vs_accuracy.png'}")

def plot_error_rates_heatmap(df: pd.DataFrame, output_dir: Path):
    """Plot error rates heatmap"""
    plt.figure(figsize=(10, 8))
    
    pivot_data = df.pivot(
        index='negative_error_rate',
        columns='positive_error_rate',
        values='accuracy'
    )
    
    sns.heatmap(
        pivot_data,
        annot=True,
        fmt='.3f',
        cmap='RdYlGn',
        vmin=0,
        vmax=1
    )
    
    plt.title('Reconstruction Accuracy vs Error Rates')
    plt.xlabel('Positive Error Rate')
    plt.ylabel('Negative Error Rate')
    plt.savefig(output_dir / 'error_rates_heatmap.png')
    plt.close()
    
    logging.info(f"Error rates heatmap saved to {output_dir / 'error_rates_heatmap.png'}")

def export_summary(length_results: pd.DataFrame, error_results: pd.DataFrame, output_dir: Path):
    """Export a summary of the benchmark results"""
    summary_file = output_dir / 'summary.txt'
    
    with open(summary_file, 'w') as f:
        # Length vs Accuracy Summary
        f.write("=== Sequence Length vs Accuracy Summary ===\n\n")
        for k in length_results['k_mer_length'].unique():
            data = length_results[length_results['k_mer_length'] == k]
            f.write(f"K-mer length = {k}:\n")
            f.write(f"  Best accuracy: {data['accuracy'].max():.3f} (length={data.loc[data['accuracy'].idxmax(), 'sequence_length']})\n")
            f.write(f"  Average accuracy: {data['accuracy'].mean():.3f} ± {data['accuracy'].std():.3f}\n")
            f.write(f"  Average execution time: {data['time'].mean():.3f}s ± {data['time'].std():.3f}s\n\n")
        
        # Error Rates Summary
        f.write("\n=== Error Rates Summary ===\n\n")
        f.write(f"Best accuracy: {error_results['accuracy'].max():.3f}\n")
        best_idx = error_results['accuracy'].idxmax()
        f.write(f"  at neg_error={error_results.loc[best_idx, 'negative_error_rate']:.2f}, "
                f"pos_error={error_results.loc[best_idx, 'positive_error_rate']:.2f}\n")
        f.write(f"Worst accuracy: {error_results['accuracy'].min():.3f}\n")
        worst_idx = error_results['accuracy'].idxmin()
        f.write(f"  at neg_error={error_results.loc[worst_idx, 'negative_error_rate']:.2f}, "
                f"pos_error={error_results.loc[worst_idx, 'positive_error_rate']:.2f}\n")
        f.write(f"Average execution time: {error_results['time'].mean():.3f}s ± {error_results['time'].std():.3f}s\n")
    
    logging.info(f"Summary exported to {summary_file}")

def main():
    # Create output directory
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir = Path('benchmark_results') / timestamp
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Setup logging
    setup_logging(output_dir)
    logging.info(f"Starting benchmarks. Results will be saved to: {output_dir}")
    
    # Run benchmarks
    logging.info("Running length vs accuracy benchmark...")
    length_results = run_length_vs_accuracy_benchmark()
    length_results.to_csv(output_dir / 'length_vs_accuracy.csv', index=False)
    plot_length_vs_accuracy(length_results, output_dir)
    
    logging.info("\nRunning error rate benchmark...")
    error_results = run_error_rate_benchmark()
    error_results.to_csv(output_dir / 'error_rates.csv', index=False)
    plot_error_rates_heatmap(error_results, output_dir)
    
    # Export summary
    export_summary(length_results, error_results, output_dir)
    
    logging.info(f"\nBenchmark results saved to: {output_dir}")

if __name__ == '__main__':
    main() 