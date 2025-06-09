import random
from typing import List, Set
from collections import Counter
from src.generators.dna_generator import DNAGenerator

class SpectrumGenerator:
    """Generator for DNA hybridization spectrum"""
    
    def __init__(self):
        self.dna_generator = DNAGenerator()
    
    def generate(self, dna: str, k: int, negative_error_rate: float = 0.0,
                positive_error_rate: float = 0.0) -> List[str]:
        """
        Generate hybridization spectrum from DNA sequence
        
        Args:
            dna: The DNA sequence
            k: Length of oligonucleotides
            negative_error_rate: Rate of k-mers to remove (false negatives)
            positive_error_rate: Rate of random k-mers to add (false positives)
            
        Returns:
            List of k-mers representing the hybridization spectrum
            
        Raises:
            ValueError: If parameters are invalid
        """
        if k < 1 or k > len(dna):
            raise ValueError("Invalid k-mer length")
            
        if not (0 <= negative_error_rate <= 1 and 0 <= positive_error_rate <= 1):
            raise ValueError("Error rates must be between 0 and 1")
        
        # Generate complete spectrum with duplicates preserved
        spectrum = []
        for i in range(len(dna) - k + 1):
            kmer = dna[i:i+k]
            spectrum.append(kmer)
        
        # Apply negative errors (remove k-mers)
        if negative_error_rate > 0:
            num_to_remove = int(len(spectrum) * negative_error_rate)
            for _ in range(num_to_remove):
                if spectrum:  # Check if there are still k-mers to remove
                    idx = random.randrange(len(spectrum))
                    spectrum.pop(idx)
        
        # Apply positive errors (add random k-mers)
        if positive_error_rate > 0:
            num_to_add = int(len(spectrum) * positive_error_rate)
            for _ in range(num_to_add):
                # Generate random k-mer
                random_kmer = self.dna_generator.generate(k)
                spectrum.append(random_kmer)
        
        return spectrum 