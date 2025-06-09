import random
from typing import List

class DNAGenerator:
    """Generator for random DNA sequences"""
    
    def __init__(self):
        self.nucleotides = ['A', 'C', 'G', 'T']
    
    def generate(self, length: int) -> str:
        """
        Generate a random DNA sequence of specified length
        
        Args:
            length: The length of DNA sequence to generate
            
        Returns:
            A string representing the DNA sequence
            
        Raises:
            ValueError: If length is less than 1
        """
        if length < 1:
            raise ValueError("DNA length must be positive")
            
        return ''.join(random.choice(self.nucleotides) for _ in range(length)) 