import pytest
from src.generators.spectrum_generator import SpectrumGenerator
from src.generators.dna_generator import DNAGenerator

def test_spectrum_generation():
    """Test basic spectrum generation without errors"""
    dna = "ACGTACGT"
    k = 3
    generator = SpectrumGenerator()
    spectrum = generator.generate(dna, k)
    
    # Expected k-mers for ACGTACGT with k=3
    expected = {'ACG', 'CGT', 'GTA', 'TAC', 'ACG', 'CGT'}
    assert set(spectrum) == expected
    assert len(spectrum) == len(dna) - k + 1

def test_spectrum_with_negative_errors():
    """Test spectrum generation with negative errors"""
    dna = "ACGTACGT"
    k = 3
    error_rate = 0.2  # 20% negative errors
    generator = SpectrumGenerator()
    spectrum = generator.generate(dna, k, negative_error_rate=error_rate)
    
    # Should have fewer k-mers than complete spectrum
    assert len(spectrum) < len(dna) - k + 1

def test_spectrum_with_positive_errors():
    """Test spectrum generation with positive errors"""
    dna = "ACGTACGT"
    k = 3
    error_rate = 0.2  # 20% positive errors
    generator = SpectrumGenerator()
    spectrum = generator.generate(dna, k, positive_error_rate=error_rate)
    
    # Should have more k-mers than complete spectrum
    assert len(spectrum) > len(dna) - k + 1

def test_invalid_k_value():
    """Test if generator raises error for invalid k value"""
    dna = "ACGTACGT"
    generator = SpectrumGenerator()
    
    with pytest.raises(ValueError):
        generator.generate(dna, 0)
    with pytest.raises(ValueError):
        generator.generate(dna, -1)
    with pytest.raises(ValueError):
        generator.generate(dna, len(dna) + 1)

def test_invalid_error_rates():
    """Test if generator raises error for invalid error rates"""
    dna = "ACGTACGT"
    k = 3
    generator = SpectrumGenerator()
    
    with pytest.raises(ValueError):
        generator.generate(dna, k, negative_error_rate=-0.1)
    with pytest.raises(ValueError):
        generator.generate(dna, k, positive_error_rate=-0.1)
    with pytest.raises(ValueError):
        generator.generate(dna, k, negative_error_rate=1.1)
    with pytest.raises(ValueError):
        generator.generate(dna, k, positive_error_rate=1.1) 