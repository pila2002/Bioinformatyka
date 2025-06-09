import time
import numpy as np
from typing import Dict, List, Tuple
from dataclasses import dataclass
from ..generators.dna_generator import DNAGenerator
from ..generators.spectrum_generator import SpectrumGenerator
from ..algorithms.classic_sbh import ClassicSBH

@dataclass
class BenchmarkParameters:
    """Parameters for SBH benchmarking"""
    sequence_length: int
    k_mer_length: int
    negative_error_rate: float
    positive_error_rate: float
    num_trials: int

@dataclass
class BenchmarkResult:
    """Results from a single benchmark run"""
    parameters: BenchmarkParameters
    avg_reconstruction_accuracy: float
    avg_execution_time: float
    std_reconstruction_accuracy: float
    std_execution_time: float
    successful_reconstructions: int

class SBHBenchmark:
    """Benchmark class for SBH algorithm"""
    
    def __init__(self):
        self.dna_generator = DNAGenerator()
        self.spectrum_generator = SpectrumGenerator()
        self.sbh = ClassicSBH()
    
    def run_single_trial(self, params: BenchmarkParameters) -> Tuple[float, float]:
        """Run a single trial with given parameters"""
        # Generate random DNA sequence
        original_dna = self.dna_generator.generate(params.sequence_length)
        
        # Generate spectrum with errors
        spectrum = self.spectrum_generator.generate(
            original_dna,
            params.k_mer_length,
            params.negative_error_rate,
            params.positive_error_rate
        )
        
        # Measure reconstruction time
        start_time = time.time()
        reconstructed_dna = self.sbh.reconstruct(
            spectrum,
            params.sequence_length,
            params.k_mer_length
        )
        execution_time = time.time() - start_time
        
        # Calculate reconstruction accuracy
        accuracy = self._calculate_accuracy(original_dna, reconstructed_dna)
        
        return accuracy, execution_time
    
    def run_benchmark(self, params: BenchmarkParameters) -> BenchmarkResult:
        """Run multiple trials and aggregate results"""
        accuracies = []
        times = []
        successful = 0
        
        for _ in range(params.num_trials):
            try:
                accuracy, execution_time = self.run_single_trial(params)
                accuracies.append(accuracy)
                times.append(execution_time)
                if accuracy == 1.0:
                    successful += 1
            except Exception as e:
                print(f"Trial failed: {e}")
                continue
        
        if not accuracies:
            raise ValueError("All trials failed")
        
        return BenchmarkResult(
            parameters=params,
            avg_reconstruction_accuracy=np.mean(accuracies),
            avg_execution_time=np.mean(times),
            std_reconstruction_accuracy=np.std(accuracies),
            std_execution_time=np.std(times),
            successful_reconstructions=successful
        )
    
    def _calculate_accuracy(self, original: str, reconstructed: str) -> float:
        """Calculate reconstruction accuracy"""
        if len(original) != len(reconstructed):
            return 0.0
        
        matches = sum(1 for a, b in zip(original, reconstructed) if a == b)
        return matches / len(original) 