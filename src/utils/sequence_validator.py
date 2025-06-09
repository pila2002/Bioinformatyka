class SequenceValidator:
    def validate_reconstruction(self, reconstructed_sequence, original_spectrum, k):
        """Validate reconstructed sequence against original spectrum."""
        # Generate spectrum from reconstructed sequence
        reconstructed_spectrum = set()
        for i in range(len(reconstructed_sequence) - k + 1):
            reconstructed_spectrum.add(reconstructed_sequence[i:i+k])
        
        # Calculate coverage
        correct_k_mers = len(reconstructed_spectrum.intersection(original_spectrum))
        total_k_mers = len(original_spectrum)
        coverage = (correct_k_mers / total_k_mers) * 100 if total_k_mers > 0 else 0
        
        # Check if all k-mers from original spectrum are present
        is_valid = reconstructed_spectrum.issuperset(original_spectrum)
        
        return is_valid, coverage 