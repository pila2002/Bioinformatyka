import pytest
from src.algorithms.classic_sbh import ClassicSBH
from src.generators.dna_generator import DNAGenerator
from src.generators.spectrum_generator import SpectrumGenerator

def test_perfect_reconstruction():
    """Test DNA reconstruction with perfect spectrum (no errors)"""
    original_dna = "ACGTACGT"
    k = 3
    spectrum_gen = SpectrumGenerator()
    spectrum = spectrum_gen.generate(original_dna, k)
    
    print(f"\nOriginal DNA: {original_dna}")
    print(f"k: {k}")
    print(f"Spectrum: {sorted(spectrum)}")
    
    sbh = ClassicSBH()
    reconstructed_dna = sbh.reconstruct(spectrum, len(original_dna), k)
    
    print(f"Reconstructed DNA: {reconstructed_dna}")
    
    assert reconstructed_dna == original_dna

def test_reconstruction_with_negative_errors():
    """Test DNA reconstruction with negative errors"""
    original_dna = "ACGTACGT"
    k = 3
    spectrum_gen = SpectrumGenerator()
    spectrum = spectrum_gen.generate(original_dna, k, negative_error_rate=0.2)
    
    sbh = ClassicSBH()
    reconstructed_dna = sbh.reconstruct(spectrum, len(original_dna), k)
    
    # With errors, we allow some length variation and check if nucleotides are valid
    assert len(reconstructed_dna) >= len(original_dna) * 0.5  # At least 50% of original length
    assert all(nucleotide in 'ACGT' for nucleotide in reconstructed_dna)
    assert len(reconstructed_dna) > 0

def test_reconstruction_with_positive_errors():
    """Test DNA reconstruction with positive errors"""
    original_dna = "ACGTACGT"
    k = 3
    spectrum_gen = SpectrumGenerator()
    spectrum = spectrum_gen.generate(original_dna, k, positive_error_rate=0.2)
    
    sbh = ClassicSBH()
    reconstructed_dna = sbh.reconstruct(spectrum, len(original_dna), k)
    
    # With errors, we allow some length variation and check if nucleotides are valid  
    assert len(reconstructed_dna) >= len(original_dna) * 0.5  # At least 50% of original length
    assert all(nucleotide in 'ACGT' for nucleotide in reconstructed_dna)
    assert len(reconstructed_dna) > 0

def test_invalid_spectrum():
    """Test reconstruction with invalid spectrum"""
    sbh = ClassicSBH()
    with pytest.raises(ValueError):
        sbh.reconstruct([], 10, 3)  # Empty spectrum
    with pytest.raises(ValueError):
        sbh.reconstruct(['ACGT'], 10, 5)  # k-mer length mismatch

def test_invalid_parameters():
    """Test reconstruction with invalid parameters"""
    spectrum = ['ACG', 'CGT']
    sbh = ClassicSBH()
    
    with pytest.raises(ValueError):
        sbh.reconstruct(spectrum, 0, 3)  # Invalid length
    with pytest.raises(ValueError):
        sbh.reconstruct(spectrum, 5, 0)  # Invalid k
    with pytest.raises(ValueError):
        sbh.reconstruct(spectrum, 2, 3)  # Length shorter than k

def test_longer_sequence():
    """Test reconstruction of a longer sequence"""
    gen = DNAGenerator()
    original_dna = gen.generate(300)
    k = 8
    
    spectrum_gen = SpectrumGenerator()
    spectrum = spectrum_gen.generate(original_dna, k)
    
    sbh = ClassicSBH()
    reconstructed_dna = sbh.reconstruct(spectrum, len(original_dna), k)
    
    assert len(reconstructed_dna) == len(original_dna)
    assert all(n in 'ACGT' for n in reconstructed_dna) 