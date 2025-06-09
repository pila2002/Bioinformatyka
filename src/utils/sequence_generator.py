import random

class SequenceGenerator:
    def __init__(self):
        self.nucleotides = ['A', 'C', 'G', 'T']
    
    def generate_random_sequence(self, length):
        """Generate a random DNA sequence of given length."""
        return ''.join(random.choice(self.nucleotides) for _ in range(length)) 