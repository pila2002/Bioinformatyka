#!/usr/bin/env python3
"""
Test script for Three-Phase SBH algorithm with adaptive strategy
Usage: python test_two_phase.py --length 300 --k 8 --error 0.05 --error_threshold 0.15 --trials 3
"""

import argparse
import time
import sys
from pathlib import Path

# Add src to path
sys.path.append(str(Path(__file__).parent / "src"))

from src.generators.dna_generator import DNAGenerator
from src.generators.spectrum_generator import SpectrumGenerator
from src.algorithms.classic_sbh import ClassicSBH
from src.algorithms.two_phase_sbh import ThreePhaseSBH

def calculate_accuracy(original: str, reconstructed: str) -> float:
    """Calculate reconstruction accuracy as percentage of correct characters"""
    if not original or not reconstructed:
        return 0.0
    
    # Handle different lengths
    min_len = min(len(original), len(reconstructed))
    if min_len == 0:
        return 0.0
    
    # Count matching characters
    matches = sum(1 for i in range(min_len) if original[i] == reconstructed[i])
    
    # Calculate accuracy based on original length
    accuracy = matches / len(original)
    return accuracy * 100

def run_comparison_test(length: int, k: int, error_rate: float, error_threshold: float, trials: int):
    """Run comparison test between Classic and Three-Phase SBH"""
    
    print(f"\n{'='*60}")
    print(f"TEST PARAMETERS:")
    print(f"Sequence length: {length}")
    print(f"K-mer size: {k}")
    print(f"Error rate: {error_rate:.1%}")
    print(f"Error threshold: {error_threshold}")
    print(f"Number of trials: {trials}")
    print(f"{'='*60}")
    
    # Initialize algorithms
    classic_sbh = ClassicSBH()
    three_phase_sbh = ThreePhaseSBH(error_threshold=error_threshold)
    
    # Results storage
    classic_results = []
    three_phase_results = []
    classic_times = []
    three_phase_times = []
    
    for trial in range(trials):
        print(f"\n--- Trial {trial + 1}/{trials} ---")
        
        # Generate test data
        dna_gen = DNAGenerator()
        original_sequence = dna_gen.generate(length)
        
        spectrum_gen = SpectrumGenerator()
        # Generate spectrum with errors directly
        spectrum = spectrum_gen.generate(original_sequence, k, error_rate, error_rate)
        
        print(f"Generated sequence of length {len(original_sequence)}")
        print(f"Spectrum size: {len(spectrum)} k-mers")
        
        # Test Classic SBH
        print(f"\nTesting Classic SBH...")
        start_time = time.time()
        classic_result = classic_sbh.reconstruct(spectrum, length, k)
        classic_time = time.time() - start_time
        classic_accuracy = calculate_accuracy(original_sequence, classic_result)
        
        classic_results.append(classic_accuracy)
        classic_times.append(classic_time)
        
        print(f"Classic SBH - Length: {len(classic_result)}, Accuracy: {classic_accuracy:.1f}%, Time: {classic_time:.3f}s")
        
        # Test Three-Phase SBH
        print(f"\nTesting Three-Phase SBH...")
        start_time = time.time()
        three_phase_result = three_phase_sbh.reconstruct(spectrum, length, k)
        three_phase_time = time.time() - start_time
        three_phase_accuracy = calculate_accuracy(original_sequence, three_phase_result)
        
        three_phase_results.append(three_phase_accuracy)
        three_phase_times.append(three_phase_time)
        
        print(f"Three-Phase SBH - Length: {len(three_phase_result)}, Accuracy: {three_phase_accuracy:.1f}%, Time: {three_phase_time:.3f}s")
        
        # Show improvement
        accuracy_improvement = three_phase_accuracy - classic_accuracy
        time_ratio = classic_time / three_phase_time if three_phase_time > 0 else float('inf')
        
        print(f"\nTrial {trial + 1} Results:")
        print(f"  Accuracy improvement: {accuracy_improvement:+.1f}%")
        print(f"  Speed ratio: {time_ratio:.1f}x {'faster' if time_ratio > 1 else 'slower'}")
    
    # Overall results
    print(f"\n{'='*60}")
    print(f"OVERALL RESULTS ({trials} trials)")
    print(f"{'='*60}")
    
    classic_avg_acc = sum(classic_results) / len(classic_results)
    three_phase_avg_acc = sum(three_phase_results) / len(three_phase_results)
    classic_avg_time = sum(classic_times) / len(classic_times)
    three_phase_avg_time = sum(three_phase_times) / len(three_phase_times)
    
    print(f"\nClassic SBH:")
    print(f"  Average accuracy: {classic_avg_acc:.1f}%")
    print(f"  Average time: {classic_avg_time:.3f}s")
    
    print(f"\nThree-Phase SBH:")
    print(f"  Average accuracy: {three_phase_avg_acc:.1f}%")
    print(f"  Average time: {three_phase_avg_time:.3f}s")
    
    print(f"\nImprovements:")
    accuracy_improvement = three_phase_avg_acc - classic_avg_acc
    time_ratio = classic_avg_time / three_phase_avg_time if three_phase_avg_time > 0 else float('inf')
    
    print(f"  Accuracy: {accuracy_improvement:+.1f}% improvement")
    print(f"  Speed: {time_ratio:.1f}x {'faster' if time_ratio > 1 else 'slower'}")
    
    # Success rate (>50% accuracy)
    classic_success_rate = len([x for x in classic_results if x > 50]) / len(classic_results) * 100
    three_phase_success_rate = len([x for x in three_phase_results if x > 50]) / len(three_phase_results) * 100
    
    print(f"\nSuccess rate (>50% accuracy):")
    print(f"  Classic SBH: {classic_success_rate:.0f}%")
    print(f"  Three-Phase SBH: {three_phase_success_rate:.0f}%")

def main():
    parser = argparse.ArgumentParser(description='Test Three-Phase SBH Algorithm')
    parser.add_argument('--length', type=int, default=200, help='DNA sequence length')
    parser.add_argument('--k', type=int, default=8, help='K-mer size')  
    parser.add_argument('--error', type=float, default=0.05, help='Error rate (0.0-1.0)')
    parser.add_argument('--error_threshold', type=float, default=0.15, help='Error threshold for adaptive strategy')
    parser.add_argument('--trials', type=int, default=1, help='Number of test trials')
    
    args = parser.parse_args()
    
    # Validate parameters
    if args.length < args.k:
        print(f"Error: Sequence length ({args.length}) must be >= k-mer size ({args.k})")
        return
    
    if not 0 <= args.error <= 1:
        print(f"Error: Error rate must be between 0.0 and 1.0")
        return
    
    if not 0 <= args.error_threshold <= 1:
        print(f"Error: Error threshold must be between 0.0 and 1.0")
        return
    
    if args.trials < 1:
        print(f"Error: Number of trials must be >= 1")
        return
    
    run_comparison_test(args.length, args.k, args.error, args.error_threshold, args.trials)

if __name__ == "__main__":
    main() 