from ..algorithms.improved_sbh import ImprovedSBH

class ImprovedSBHBenchmark(SBHBenchmark):
    """Benchmark for improved SBH algorithm"""
    
    def __init__(self):
        super().__init__()
        self.sbh = ImprovedSBH()  # Use improved algorithm
        
    # Rest of the methods remain the same 