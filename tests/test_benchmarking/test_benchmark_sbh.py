import pytest
from src.benchmarking.benchmark_sbh import SBHBenchmark, BenchmarkParameters

def test_benchmark_perfect_reconstruction():
    """Test benchmarking with perfect reconstruction"""
    benchmark = SBHBenchmark()
    params = BenchmarkParameters(
        sequence_length=100,
        k_mer_length=8,
        negative_error_rate=0.0,
        positive_error_rate=0.0,
        num_trials=3
    )
    
    result = benchmark.run_benchmark(params)
    
    assert result.avg_reconstruction_accuracy > 0.9  # Should be near perfect
    assert result.avg_execution_time > 0
    assert result.std_reconstruction_accuracy >= 0
    assert result.std_execution_time >= 0
    assert result.successful_reconstructions > 0

def test_benchmark_with_errors():
    """Test benchmarking with error rates"""
    benchmark = SBHBenchmark()
    params = BenchmarkParameters(
        sequence_length=100,
        k_mer_length=8,
        negative_error_rate=0.1,
        positive_error_rate=0.1,
        num_trials=3
    )
    
    result = benchmark.run_benchmark(params)
    
    assert 0 <= result.avg_reconstruction_accuracy <= 1
    assert result.avg_execution_time > 0
    assert result.std_reconstruction_accuracy >= 0
    assert result.std_execution_time >= 0

def test_invalid_parameters():
    """Test benchmarking with invalid parameters"""
    benchmark = SBHBenchmark()
    
    # Invalid sequence length
    with pytest.raises(ValueError):
        params = BenchmarkParameters(
            sequence_length=0,
            k_mer_length=8,
            negative_error_rate=0.0,
            positive_error_rate=0.0,
            num_trials=3
        )
        benchmark.run_benchmark(params)
    
    # Invalid k-mer length
    with pytest.raises(ValueError):
        params = BenchmarkParameters(
            sequence_length=100,
            k_mer_length=0,
            negative_error_rate=0.0,
            positive_error_rate=0.0,
            num_trials=3
        )
        benchmark.run_benchmark(params)
    
    # Invalid error rate
    with pytest.raises(ValueError):
        params = BenchmarkParameters(
            sequence_length=100,
            k_mer_length=8,
            negative_error_rate=1.5,  # > 1
            positive_error_rate=0.0,
            num_trials=3
        )
        benchmark.run_benchmark(params)

def test_accuracy_calculation():
    """Test accuracy calculation"""
    benchmark = SBHBenchmark()
    
    # Perfect match
    assert benchmark._calculate_accuracy("ACGT", "ACGT") == 1.0
    
    # No match
    assert benchmark._calculate_accuracy("AAAA", "TTTT") == 0.0
    
    # Partial match
    assert benchmark._calculate_accuracy("ACGT", "ACTT") == 0.75
    
    # Different lengths
    assert benchmark._calculate_accuracy("ACGT", "ACG") == 0.0 