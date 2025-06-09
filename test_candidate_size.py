#!/usr/bin/env python3
"""
Test script for candidate_size parameter with Levenshtein similarity
Tests with fixed parameters: n=300, k=8, 5% positive and negative errors
Candidate sizes: {5, 10, 15, 20, 25} with 3 repetitions each
"""

import sys
import time
import csv
from pathlib import Path
from typing import List, Tuple

# Add src to path
sys.path.append(str(Path(__file__).parent / "src"))

from src.generators.dna_generator import DNAGenerator
from src.generators.spectrum_generator import SpectrumGenerator
from src.algorithms.two_phase_sbh import ThreePhaseSBH

def levenshtein_distance(s1: str, s2: str) -> int:
    """Calculate Levenshtein distance between two strings"""
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    
    return previous_row[-1]

def levenshtein_similarity(s1: str, s2: str) -> float:
    """Calculate Levenshtein similarity as percentage"""
    if not s1 and not s2:
        return 100.0
    
    distance = levenshtein_distance(s1, s2)
    max_length = max(len(s1), len(s2))
    
    if max_length == 0:
        return 100.0
    
    similarity = (1 - distance / max_length) * 100
    return max(0.0, similarity)  # Ensure non-negative

def run_candidate_size_test():
    """Run systematic test of candidate_size parameter"""
    
    # Test parameters
    SEQUENCE_LENGTH = 300
    K_MER_SIZE = 8
    POSITIVE_ERROR_RATE = 0.05
    NEGATIVE_ERROR_RATE = 0.05
    CANDIDATE_SIZES = [5, 10, 15, 20, 25]
    REPETITIONS = 3
    
    print("="*80)
    print("CANDIDATE SIZE TEST WITH LEVENSHTEIN SIMILARITY")
    print("="*80)
    print(f"Sequence length: {SEQUENCE_LENGTH}")
    print(f"K-mer size: {K_MER_SIZE}")
    print(f"Positive error rate: {POSITIVE_ERROR_RATE:.1%}")
    print(f"Negative error rate: {NEGATIVE_ERROR_RATE:.1%}")
    print(f"Candidate sizes: {CANDIDATE_SIZES}")
    print(f"Repetitions per size: {REPETITIONS}")
    print("="*80)
    
    # Results storage
    results = []
    
    # Initialize generators
    dna_gen = DNAGenerator()
    spectrum_gen = SpectrumGenerator()
    
    for candidate_size in CANDIDATE_SIZES:
        print(f"\n--- Testing candidate_size = {candidate_size} ---")
        
        candidate_results = []
        
        for repetition in range(REPETITIONS):
            print(f"  Repetition {repetition + 1}/{REPETITIONS}")
            
            # Generate test data
            original_sequence = dna_gen.generate(SEQUENCE_LENGTH)
            spectrum = spectrum_gen.generate(
                original_sequence, 
                K_MER_SIZE, 
                NEGATIVE_ERROR_RATE, 
                POSITIVE_ERROR_RATE
            )
            
            print(f"    Generated sequence length: {len(original_sequence)}")
            print(f"    Spectrum size: {len(spectrum)} k-mers")
            
            # Test algorithm with current candidate_size
            algorithm = ThreePhaseSBH(
                error_threshold=0.15, 
                candidate_size=candidate_size
            )
            
            start_time = time.time()
            reconstructed = algorithm.reconstruct(spectrum, SEQUENCE_LENGTH, K_MER_SIZE)
            execution_time = time.time() - start_time
            
            # Calculate Levenshtein similarity
            similarity = levenshtein_similarity(original_sequence, reconstructed)
            
            # Store results
            result = {
                'candidate_size': candidate_size,
                'repetition': repetition + 1,
                'original_length': len(original_sequence),
                'reconstructed_length': len(reconstructed),
                'spectrum_size': len(spectrum),
                'levenshtein_similarity': similarity,
                'execution_time': execution_time
            }
            
            results.append(result)
            candidate_results.append(similarity)
            
            print(f"    Reconstructed length: {len(reconstructed)}")
            print(f"    Levenshtein similarity: {similarity:.2f}%")
            print(f"    Execution time: {execution_time:.3f}s")
        
        # Calculate statistics for this candidate_size
        avg_similarity = sum(candidate_results) / len(candidate_results)
        min_similarity = min(candidate_results)
        max_similarity = max(candidate_results)
        
        print(f"  Summary for candidate_size {candidate_size}:")
        print(f"    Average similarity: {avg_similarity:.2f}%")
        print(f"    Min similarity: {min_similarity:.2f}%")
        print(f"    Max similarity: {max_similarity:.2f}%")
    
    # Save results to CSV
    csv_filename = f"candidate_size_test_results_{int(time.time())}.csv"
    with open(csv_filename, 'w', newline='') as csvfile:
        fieldnames = [
            'candidate_size', 'repetition', 'original_length', 
            'reconstructed_length', 'spectrum_size', 
            'levenshtein_similarity', 'execution_time'
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    
    print(f"\n--- DETAILED RESULTS ---")
    print(f"Results saved to: {csv_filename}")
    
    # Calculate and display summary statistics
    print(f"\n--- SUMMARY STATISTICS ---")
    print(f"{'Candidate Size':<15} {'Avg Similarity':<15} {'Min Similarity':<15} {'Max Similarity':<15} {'Avg Time':<15}")
    print("-" * 75)
    
    for candidate_size in CANDIDATE_SIZES:
        size_results = [r for r in results if r['candidate_size'] == candidate_size]
        
        avg_sim = sum(r['levenshtein_similarity'] for r in size_results) / len(size_results)
        min_sim = min(r['levenshtein_similarity'] for r in size_results)
        max_sim = max(r['levenshtein_similarity'] for r in size_results)
        avg_time = sum(r['execution_time'] for r in size_results) / len(size_results)
        
        print(f"{candidate_size:<15} {avg_sim:<15.2f} {min_sim:<15.2f} {max_sim:<15.2f} {avg_time:<15.3f}")
    
    # Find best candidate_size
    best_candidate_size = None
    best_avg_similarity = -1
    
    for candidate_size in CANDIDATE_SIZES:
        size_results = [r for r in results if r['candidate_size'] == candidate_size]
        avg_sim = sum(r['levenshtein_similarity'] for r in size_results) / len(size_results)
        
        if avg_sim > best_avg_similarity:
            best_avg_similarity = avg_sim
            best_candidate_size = candidate_size
    
    print(f"\n--- CONCLUSION ---")
    print(f"Best candidate_size: {best_candidate_size}")
    print(f"Best average similarity: {best_avg_similarity:.2f}%")
    
    return results

if __name__ == "__main__":
    results = run_candidate_size_test() 