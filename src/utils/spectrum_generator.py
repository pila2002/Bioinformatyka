class SpectrumGenerator:
    def generate_spectrum(self, sequence, k):
        """Generate k-mer spectrum from a sequence."""
        spectrum = set()
        for i in range(len(sequence) - k + 1):
            spectrum.add(sequence[i:i+k])
        return spectrum 