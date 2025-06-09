from typing import List, Dict, Set, Tuple, Optional
import networkx as nx
from collections import defaultdict, Counter
import numpy as np
from dataclasses import dataclass
import time

@dataclass
class Contig:
    """Represents a contiguous sequence fragment"""
    sequence: str
    confidence: float
    start_kmer: str
    end_kmer: str
    supporting_kmers: Set[str]

class ThreePhaseSBH:
    """Three-phase adaptive SBH algorithm with rescue mechanisms"""
    
    def __init__(self, error_threshold: float = 0.15, candidate_size: int = 10):
        self._spectrum_analysis = {}
        self._contig_graph = None
        self._phase1_confidence_threshold = 0.6
        self._min_contig_length = None
        self.error_threshold = error_threshold  # Adaptive strategy threshold
        self.candidate_size = candidate_size  # Number of candidates to consider
        self._used_kmers = set()
        self._adaptive_strategy = "conservative"  # conservative, aggressive, rescue
    
    def reconstruct(self, spectrum: List[str], target_length: int, k: int) -> str:
        """
        Three-phase reconstruction with adaptive strategy:
        Phase 1: Build reliable contigs
        Phase 2: Connect contigs and extend greedily  
        Phase 3: Rescue jumps to unvisited graph regions
        """
        print(f"\n=== Three-Phase SBH Reconstruction ===")
        print(f"Target: {target_length}, k: {k}, Spectrum: {len(spectrum)} k-mers")
        print(f"Error threshold: {self.error_threshold}")
        
        self._k = k
        self._target_length = target_length
        self._used_kmers = set()
        
        # Analyze spectrum quality to determine strategy
        self._analyze_spectrum_quality(spectrum)
        
        # Phase 1: Build reliable contigs
        print(f"\n--- Phase 1: Building Contigs ---")
        contigs = self._build_contigs(spectrum, k)
        print(f"Built {len(contigs)} contigs")
        
        if contigs:
            # Phase 2: Connect contigs
            print(f"\n--- Phase 2: Connecting Contigs ---")
            sequence = self._connect_contigs(contigs, spectrum, k, target_length)
        else:
            print("No contigs built, starting with greedy approach")
            sequence = ""
        
        # Phase 3: Rescue mechanism with graph jumping
        print(f"\n--- Phase 3: Rescue Reconstruction ---")
        final_sequence = self._rescue_reconstruction(sequence, spectrum, k, target_length)
        
        return final_sequence
    
    def _analyze_spectrum_quality(self, spectrum: List[str]):
        """Analyze spectrum to determine adaptive strategy"""
        unique_kmers = len(set(spectrum))
        total_kmers = len(spectrum)
        coverage = total_kmers / unique_kmers if unique_kmers > 0 else 1
        
        # Calculate k-mer frequency distribution
        kmer_counts = Counter(spectrum)
        freq_variance = np.var(list(kmer_counts.values()))
        
        # Determine strategy based on quality metrics
        if coverage < 1.2 and freq_variance < 2:
            self._adaptive_strategy = "aggressive"  # High quality data
            self._phase1_confidence_threshold = 0.5
        elif coverage > 2.0 or freq_variance > 10:
            self._adaptive_strategy = "rescue"  # Poor quality data  
            self._phase1_confidence_threshold = 0.8
        else:
            self._adaptive_strategy = "conservative"  # Medium quality
            self._phase1_confidence_threshold = 0.6
        
        print(f"Spectrum quality analysis:")
        print(f"  Coverage: {coverage:.2f}, Variance: {freq_variance:.2f}")
        print(f"  Adaptive strategy: {self._adaptive_strategy}")
    
    def _build_contigs(self, spectrum: List[str], k: int) -> List[Contig]:
        """Phase 1: Build reliable contigs from high-confidence k-mers"""
        # Analyze k-mer reliability
        reliable_kmers = self._identify_reliable_kmers(spectrum)
        print(f"Identified {len(reliable_kmers)} reliable k-mers out of {len(spectrum)}")
        
        # Build overlap graph for reliable k-mers only
        overlap_graph = self._build_overlap_graph(reliable_kmers, k-1)
        
        # Find simple paths (contigs) in the graph
        contigs = self._extract_contigs(overlap_graph, reliable_kmers)
        
        return contigs
    
    def _identify_reliable_kmers(self, spectrum: List[str]) -> Set[str]:
        """Identify k-mers that are likely to be correct"""
        kmer_counts = Counter(spectrum)
        reliable_kmers = set()
        
        # Adjust thresholds based on adaptive strategy
        if self._adaptive_strategy == "aggressive":
            min_freq = 1
            max_freq = max(kmer_counts.values()) * 2
        elif self._adaptive_strategy == "rescue":
            min_freq = 2
            max_freq = max(kmer_counts.values()) * 0.8
        else:  # conservative
            min_freq = 1
            max_freq = max(kmer_counts.values()) * 1.5
        
        for kmer, count in kmer_counts.items():
            if min_freq <= count <= max_freq:
                # Additional quality checks
                if self._assess_kmer_quality(kmer, spectrum):
                    reliable_kmers.add(kmer)
        
        return reliable_kmers
    
    def _assess_kmer_quality(self, kmer: str, spectrum: List[str]) -> bool:
        """Assess if a k-mer is likely to be correct"""
        # Check for low-complexity regions
        if len(set(kmer)) < len(kmer) * 0.5:  # Too repetitive
            return False
        
        # Check neighborhood consistency
        prefix = kmer[:-1]
        suffix = kmer[1:]
        
        prefix_neighbors = [s for s in spectrum if s.startswith(prefix)]
        suffix_neighbors = [s for s in spectrum if s.endswith(suffix)]
        
        # Should have reasonable number of neighbors
        return len(prefix_neighbors) >= 1 and len(suffix_neighbors) >= 1
    
    def _build_overlap_graph(self, kmers: Set[str], overlap_len: int) -> nx.DiGraph:
        """Build directed overlap graph"""
        G = nx.DiGraph()
        
        for kmer in kmers:
            G.add_node(kmer)
        
        for kmer1 in kmers:
            suffix = kmer1[-overlap_len:]
            for kmer2 in kmers:
                if kmer1 != kmer2 and kmer2.startswith(suffix):
                    # Add edge with weight based on reliability
                    weight = self._calculate_edge_weight(kmer1, kmer2)
                    G.add_edge(kmer1, kmer2, weight=weight)
        
        return G
    
    def _calculate_edge_weight(self, kmer1: str, kmer2: str) -> float:
        """Calculate edge weight based on connection reliability"""
        # Simple weight based on k-mer quality
        base_weight = 1.0
        
        # Penalize low-complexity connections
        overlap = kmer1[1:]
        if len(set(overlap)) < len(overlap) * 0.6:
            base_weight *= 0.5
        
        return base_weight
    
    def _extract_contigs(self, graph: nx.DiGraph, kmers: Set[str]) -> List[Contig]:
        """Extract contigs from overlap graph"""
        contigs = []
        visited = set()
        
        # Find simple paths (no branching)
        for start_node in graph.nodes():
            if start_node in visited:
                continue
                
            if graph.in_degree(start_node) == 0 or graph.in_degree(start_node) > 1:
                # Potential start of contig
                path = self._follow_simple_path(graph, start_node, visited)
                if len(path) >= 2:  # At least 2 k-mers
                    contig_seq = self._path_to_sequence(path)
                    confidence = self._calculate_contig_confidence(path)
                    
                    contig = Contig(
                        sequence=contig_seq,
                        confidence=confidence,
                        start_kmer=path[0],
                        end_kmer=path[-1],
                        supporting_kmers=set(path)
                    )
                    contigs.append(contig)
        
        return contigs
    
    def _follow_simple_path(self, graph: nx.DiGraph, start: str, visited: set) -> List[str]:
        """Follow a simple path without branching"""
        path = [start]
        visited.add(start)
        current = start
        
        while True:
            successors = list(graph.successors(current))
            # Continue only if exactly one unvisited successor
            unvisited_successors = [s for s in successors if s not in visited]
            
            if len(unvisited_successors) == 1:
                next_node = unvisited_successors[0]
                # Check if next node has only one predecessor (no convergence)
                if graph.in_degree(next_node) == 1:
                    path.append(next_node)
                    visited.add(next_node)
                    current = next_node
                else:
                    break
            else:
                break
        
        return path
    
    def _path_to_sequence(self, path: List[str]) -> str:
        """Convert k-mer path to sequence"""
        if not path:
            return ""
        
        sequence = path[0]
        for kmer in path[1:]:
            sequence += kmer[-1]  # Add last character
        
        return sequence
    
    def _calculate_contig_confidence(self, path: List[str]) -> float:
        """Calculate confidence score for a contig"""
        if not path:
            return 0.0
        
        # Base confidence on path length and k-mer quality
        length_bonus = min(len(path) / 10, 1.0)  # Longer contigs are better
        base_confidence = 0.5 + length_bonus * 0.3
        
        return base_confidence
    
    def _connect_contigs(self, contigs: List[Contig], spectrum: List[str], k: int, target_length: int) -> str:
        """Phase 2: Connect contigs and extend greedily"""
        if not contigs:
            return ""
        
        # Start with highest confidence contig
        contigs.sort(key=lambda c: c.confidence, reverse=True)
        current_sequence = contigs[0].sequence
        self._used_kmers.update(contigs[0].supporting_kmers)
        
        # Try to extend with other contigs
        remaining_contigs = contigs[1:]
        
        while len(current_sequence) < target_length and remaining_contigs:
            extended = False
            
            for i, contig in enumerate(remaining_contigs):
                # Try to connect at the end
                connection = self._find_connection(current_sequence, contig.sequence, k)
                if connection:
                    current_sequence += connection
                    self._used_kmers.update(contig.supporting_kmers)
                    remaining_contigs.pop(i)
                    extended = True
                    break
            
            if not extended:
                break
        
        return current_sequence
    
    def _find_connection(self, seq1: str, seq2: str, k: int) -> Optional[str]:
        """Find connection between two sequences"""
        # Look for overlap
        for overlap_len in range(k-1, 0, -1):
            if seq1.endswith(seq2[:overlap_len]):
                return seq2[overlap_len:]
        
        return None
    
    def _rescue_reconstruction(self, initial_sequence: str, spectrum: List[str], k: int, target_length: int) -> str:
        """Phase 3: Rescue reconstruction with graph jumping"""
        if not initial_sequence:
            # Start with best k-mer if no initial sequence
            initial_sequence = self._find_best_starting_kmer(spectrum, k)
        
        current_sequence = initial_sequence
        available_kmers = [s for s in spectrum if s not in self._used_kmers]
        
        print(f"Starting rescue with sequence length: {len(current_sequence)}")
        print(f"Available k-mers: {len(available_kmers)}")
        
        max_iterations = target_length * 2  # Prevent infinite loops
        iteration = 0
        
        while len(current_sequence) < target_length and available_kmers and iteration < max_iterations:
            iteration += 1
            
            # Phase 3a: Try standard extension
            extended = self._try_standard_extension(current_sequence, available_kmers, k)
            
            if extended:
                current_sequence = extended
                # Remove used k-mers
                available_kmers = [s for s in available_kmers if s not in self._get_sequence_kmers(current_sequence, k)]
                continue
            
            # Phase 3b: Try adaptive jump
            jump_result = self._try_adaptive_jump(current_sequence, available_kmers, k, target_length)
            
            if jump_result:
                current_sequence = jump_result
                available_kmers = [s for s in available_kmers if s not in self._get_sequence_kmers(current_sequence, k)]
                continue
            
            # Phase 3c: Emergency connection
            emergency_result = self._emergency_connection(current_sequence, available_kmers, k)
            
            if emergency_result:
                current_sequence = emergency_result
                available_kmers = [s for s in available_kmers if s not in self._get_sequence_kmers(current_sequence, k)]
                continue
            
            # If nothing works, break to avoid infinite loop
            break
        
        print(f"Rescue completed after {iteration} iterations")
        print(f"Final sequence length: {len(current_sequence)}/{target_length}")
        
        return current_sequence[:target_length] if len(current_sequence) > target_length else current_sequence
    
    def _find_best_starting_kmer(self, spectrum: List[str], k: int) -> str:
        """Find the best k-mer to start reconstruction"""
        kmer_scores = {}
        
        for kmer in set(spectrum):
            score = self._calculate_starting_score(kmer, spectrum, k)
            kmer_scores[kmer] = score
        
        best_kmer = max(kmer_scores.keys(), key=lambda x: kmer_scores[x])
        self._used_kmers.add(best_kmer)
        
        return best_kmer
    
    def _calculate_starting_score(self, kmer: str, spectrum: List[str], k: int) -> float:
        """Calculate score for starting k-mer"""
        # Prefer k-mers with good connectivity
        prefix = kmer[:-1]
        suffix = kmer[1:]
        
        out_degree = len([s for s in spectrum if s.startswith(suffix)])
        in_degree = len([s for s in spectrum if s.endswith(prefix)])
        
        # Prefer k-mers that are potential starts (low in-degree) or have good connectivity
        connectivity_score = out_degree - in_degree * 0.5
        complexity_score = len(set(kmer)) / len(kmer)  # Prefer complex k-mers
        
        return connectivity_score + complexity_score
    
    def _try_standard_extension(self, sequence: str, available_kmers: List[str], k: int) -> Optional[str]:
        """Try standard greedy extension with candidate_size limit"""
        if len(sequence) < k:
            return None
        
        current_suffix = sequence[-(k-1):]
        candidates = []
        
        # Find all matching k-mers
        for kmer in available_kmers:
            if kmer.startswith(current_suffix):
                score = self._calculate_extension_score(kmer, available_kmers, k)
                candidates.append((score, kmer))
        
        if not candidates:
            return None
        
        # Sort by score and take top candidate_size candidates
        candidates.sort(reverse=True)
        top_candidates = candidates[:self.candidate_size]
        
        # Choose the best one (first in sorted list)
        best_score, best_kmer = top_candidates[0]
        
        if best_kmer:
            return sequence + best_kmer[-1]
        
        return None
    
    def _calculate_extension_score(self, kmer: str, available_kmers: List[str], k: int) -> float:
        """Calculate score for extending with a k-mer"""
        # Base score
        score = 1.0
        
        # Bonus for k-mers that can be further extended
        suffix = kmer[1:]
        next_options = len([s for s in available_kmers if s.startswith(suffix) and s != kmer])
        score += next_options * 0.1
        
        # Bonus for complexity
        complexity = len(set(kmer)) / len(kmer)
        score += complexity
        
        return score
    
    def _try_adaptive_jump(self, sequence: str, available_kmers: List[str], k: int, target_length: int) -> Optional[str]:
        """Try adaptive jumping to unvisited graph regions"""
        if not available_kmers:
            return None
        
        # Find best jump target based on strategy
        if self._adaptive_strategy == "aggressive":
            jump_target = self._find_aggressive_jump(sequence, available_kmers, k)
        elif self._adaptive_strategy == "rescue":
            jump_target = self._find_rescue_jump(sequence, available_kmers, k)
        else:  # conservative
            jump_target = self._find_conservative_jump(sequence, available_kmers, k)
        
        if jump_target:
            # Add gap filler and jump
            gap_length = max(1, min(5, target_length - len(sequence) - len(jump_target)))
            gap = 'N' * gap_length  # Use N's to mark gaps
            
            return sequence + gap + jump_target
        
        return None
    
    def _find_aggressive_jump(self, sequence: str, available_kmers: List[str], k: int) -> Optional[str]:
        """Aggressive jumping strategy - prefer high-connectivity k-mers with candidate_size limit"""
        candidates = []
        
        for kmer in available_kmers:
            # Score based on connectivity potential
            suffix = kmer[1:]
            connectivity = len([s for s in available_kmers if s.startswith(suffix)])
            candidates.append((connectivity, kmer))
        
        if not candidates:
            return None
        
        # Sort by connectivity and take top candidates
        candidates.sort(reverse=True)
        top_candidates = candidates[:self.candidate_size]
        
        # Return the best one
        return top_candidates[0][1] if top_candidates else None
    
    def _find_conservative_jump(self, sequence: str, available_kmers: List[str], k: int) -> Optional[str]:
        """Conservative jumping - look for partial overlaps"""
        current_end = sequence[-(k-2):] if len(sequence) >= k-2 else sequence
        
        best_kmer = None
        best_overlap = -1
        
        for kmer in available_kmers:
            # Look for partial overlaps
            for i in range(min(len(current_end), k-1), 0, -1):
                if current_end.endswith(kmer[:i]):
                    if i > best_overlap:
                        best_overlap = i
                        best_kmer = kmer
                    break
        
        return best_kmer
    
    def _find_rescue_jump(self, sequence: str, available_kmers: List[str], k: int) -> Optional[str]:
        """Rescue jumping - just pick any reasonable k-mer"""
        if not available_kmers:
            return None
        
        # Filter out low-complexity k-mers
        good_kmers = [kmer for kmer in available_kmers 
                     if len(set(kmer)) >= len(kmer) * 0.6]
        
        if good_kmers:
            # Pick one with good extension potential
            return max(good_kmers, key=lambda x: len([s for s in available_kmers if s.startswith(x[1:])]))
        else:
            return available_kmers[0]  # Last resort
    
    def _emergency_connection(self, sequence: str, available_kmers: List[str], k: int) -> Optional[str]:
        """Emergency connection when all else fails"""
        if not available_kmers:
            return None
        
        # Just pick the first available k-mer and connect with gap
        emergency_kmer = available_kmers[0]
        gap = 'N' * min(3, k-1)  # Small gap
        
        return sequence + gap + emergency_kmer
    
    def _get_sequence_kmers(self, sequence: str, k: int) -> Set[str]:
        """Get all k-mers in a sequence"""
        if len(sequence) < k:
            return set()
        
        return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}

# Keep backward compatibility
TwoPhaseSBH = ThreePhaseSBH 