import pytest
from src.generators.dna_generator import DNAGenerator

def test_dna_generator_length():
    """Test if generated DNA has correct length"""
    length = 300
    generator = DNAGenerator()
    dna = generator.generate(length)
    assert len(dna) == length

def test_dna_generator_valid_nucleotides():
    """Test if generated DNA contains only valid nucleotides"""
    generator = DNAGenerator()
    dna = generator.generate(100)
    valid_nucleotides = set('ACGT')
    assert all(n in valid_nucleotides for n in dna)

def test_dna_generator_random():
    """Test if generator produces different sequences"""
    generator = DNAGenerator()
    dna1 = generator.generate(100)
    dna2 = generator.generate(100)
    assert dna1 != dna2  # Very unlikely to be equal if truly random

def test_dna_generator_invalid_length():
    """Test if generator raises error for invalid length"""
    generator = DNAGenerator()
    with pytest.raises(ValueError):
        generator.generate(-1)
    with pytest.raises(ValueError):
        generator.generate(0) 