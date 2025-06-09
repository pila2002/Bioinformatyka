from typing import List, Dict, Set, Tuple, Optional, Counter as CounterType
import networkx as nx
from collections import defaultdict, Counter
import numpy as np
from itertools import combinations
import heapq

class ImprovedSBH:
    """Improved Sequencing by Hybridization algorithm with enhanced accuracy"""
    
    def __init__(self):
        self._graph = None
        self._kmer_reliability = {}
        self._error_threshold = 0.3
    
    def reconstruct(self, spectrum: List[str], target_length: int, k: int) -> str:
        """
        Reconstruct DNA sequence using improved SBH algorithm
        
        Key improvements:
        1. Error detection and correction in spectrum
        2. Weighted graph construction
        3. Multiple path exploration
        4. Consensus-based reconstruction
        """
        if not spectrum or target_length < k or k < 1:
            raise ValueError("Invalid parameters")
            
        print(f"\n=== Improved SBH Reconstruction ===")
        print(f"Target: {target_length}, k: {k}, Spectrum: {len(spectrum)} k-mers")
        
        # Step 1: Preprocess spectrum to detect and correct errors
        corrected_spectrum = self._preprocess_spectrum(spectrum, k)
        print(f"Corrected spectrum: {len(corrected_spectrum)} k-mers")
        
        # Step 2: Build weighted graph with reliability scores
        self._graph = self._build_weighted_graph(corrected_spectrum, k)
        print(f"Graph: {self._graph.number_of_nodes()} nodes, {self._graph.number_of_edges()} edges")
        
        # Step 3: Find multiple candidate paths
        candidate_paths = self._find_multiple_paths(self._graph, target_length, k)
        print(f"Found {len(candidate_paths)} candidate paths")
        
        # Step 4: Score and select best path
        if candidate_paths:
            best_sequence = self._select_best_sequence(candidate_paths, corrected_spectrum, target_length, k)
            return best_sequence
        
        # Fallback to enhanced greedy reconstruction
        return self._enhanced_greedy_reconstruction(corrected_spectrum, target_length, k)
    
    def _preprocess_spectrum(self, spectrum: List[str], k: int) -> List[str]:
        """Detect and correct errors in the spectrum"""
        kmer_counts = Counter(spectrum)
        
        # 1. Remove k-mers that appear too infrequently (likely errors)
        total_kmers = len(spectrum)
        min_count = max(1, int(0.01 * total_kmers))  # At least 1% occurrence
        
        # 2. Detect hamming distance-1 k-mers that could be corrections
        reliable_kmers = {}
        questionable_kmers = {}
        
        for kmer, count in kmer_counts.items():
            if count >= min_count:
                reliable_kmers[kmer] = count
            else:
                questionable_kmers[kmer] = count
        
        # 3. Try to correct questionable k-mers
        corrected_spectrum = list(reliable_kmers.keys())
        
        for q_kmer, q_count in questionable_kmers.items():
            best_correction = None
            min_distance = float('inf')
            
            for r_kmer in reliable_kmers:
                distance = self._hamming_distance(q_kmer, r_kmer)
                if distance == 1 and distance < min_distance:
                    min_distance = distance
                    best_correction = r_kmer
            
            if best_correction:
                # Add corrected k-mer with adjusted count
                for _ in range(q_count):
                    corrected_spectrum.append(best_correction)
            else:
                # Keep original if no good correction found
                for _ in range(q_count):
                    corrected_spectrum.append(q_kmer)
        
        return corrected_spectrum
    
    def _hamming_distance(self, s1: str, s2: str) -> int:
        """Calculate Hamming distance between two strings"""
        if len(s1) != len(s2):
            return float('inf')
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))
    
    def _build_weighted_graph(self, spectrum: List[str], k: int) -> nx.DiGraph:
        """Build weighted de Bruijn graph with reliability scores"""
        kmer_counts = Counter(spectrum)
        graph = nx.DiGraph()
        
        # Calculate reliability scores for each k-mer
        total_count = sum(kmer_counts.values())
        avg_count = total_count / len(kmer_counts)
        
        for kmer, count in kmer_counts.items():
            # Reliability based on frequency and neighborhood consistency
            reliability = min(1.0, count / avg_count)
            
            # Check neighborhood consistency (how well k-mer fits with neighbors)
            prefix = kmer[:-1]
            suffix = kmer[1:]
            
            prefix_neighbors = [k for k in kmer_counts if k[:-1] == prefix]
            suffix_neighbors = [k for k in kmer_counts if k[1:] == suffix]
            
            neighborhood_score = (len(prefix_neighbors) + len(suffix_neighbors)) / 10.0
            reliability *= min(1.0, neighborhood_score)
            
            self._kmer_reliability[kmer] = reliability
            
            # Add edge with weight proportional to reliability and count
            weight = reliability * count
            graph.add_edge(prefix, suffix, weight=weight, kmer=kmer, reliability=reliability)
        
        return graph
    
    def _find_multiple_paths(self, graph: nx.DiGraph, target_length: int, k: int) -> List[List[str]]:
        """Find multiple candidate paths using different strategies"""
        paths = []
        
        # Strategy 1: Highest weight path
        try:
            path1 = self._find_maximum_weight_path(graph, target_length, k)
            if path1:
                paths.append(path1)
        except:
            pass
        
        # Strategy 2: Eulerian-based path
        try:
            path2 = self._find_modified_eulerian_path(graph, target_length, k)
            if path2:
                paths.append(path2)
        except:
            pass
        
        # Strategy 3: Multiple random walks with different starting points
        start_nodes = self._find_good_start_nodes(graph, 3)
        for start_node in start_nodes:
            try:
                path = self._guided_random_walk(graph, start_node, target_length, k)
                if path:
                    paths.append(path)
            except:
                continue
        
        return paths
    
    def _find_maximum_weight_path(self, graph: nx.DiGraph, target_length: int, k: int) -> Optional[List[str]]:
        """Find path that maximizes edge weights"""
        start_nodes = self._find_good_start_nodes(graph, 1)
        if not start_nodes:
            return None
        
        start_node = start_nodes[0]
        path = [start_node]
        used_edges = set()
        
        current_node = start_node
        while len(''.join(path)) < target_length:
            # Find best next edge
            best_edge = None
            best_weight = -1
            
            for neighbor in graph.neighbors(current_node):
                edge_data = graph[current_node][neighbor]
                edge_key = (current_node, neighbor)
                
                if edge_key not in used_edges:
                    weight = edge_data.get('weight', 0)
                    reliability = edge_data.get('reliability', 0)
                    combined_score = weight * reliability
                    
                    if combined_score > best_weight:
                        best_weight = combined_score
                        best_edge = (neighbor, edge_key)
            
            if not best_edge:
                break
                
            next_node, edge_key = best_edge
            used_edges.add(edge_key)
            path.append(next_node)
            current_node = next_node
        
        return path if len(''.join(path)) >= target_length * 0.8 else None
    
    def _find_good_start_nodes(self, graph: nx.DiGraph, max_nodes: int) -> List[str]:
        """Find good starting nodes based on degree and edge weights"""
        node_scores = {}
        
        for node in graph.nodes():
            # Score based on out-degree and outgoing edge weights
            out_degree = graph.out_degree(node)
            in_degree = graph.in_degree(node)
            
            weight_sum = sum(graph[node][neighbor].get('weight', 0) 
                           for neighbor in graph.neighbors(node))
            
            # Prefer nodes with more outgoing edges and higher weights
            score = out_degree * 2 + weight_sum - abs(out_degree - in_degree)
            node_scores[node] = score
        
        # Return top scoring nodes
        sorted_nodes = sorted(node_scores.items(), key=lambda x: x[1], reverse=True)
        return [node for node, score in sorted_nodes[:max_nodes]]
    
    def _guided_random_walk(self, graph: nx.DiGraph, start_node: str, target_length: int, k: int) -> Optional[List[str]]:
        """Perform guided random walk favoring high-weight edges"""
        path = [start_node]
        used_edges = set()
        current_node = start_node
        
        while len(''.join(path)) < target_length:
            neighbors = list(graph.neighbors(current_node))
            if not neighbors:
                break
            
            # Calculate probabilities based on edge weights
            weights = []
            valid_neighbors = []
            
            for neighbor in neighbors:
                edge_key = (current_node, neighbor)
                if edge_key not in used_edges:
                    weight = graph[current_node][neighbor].get('weight', 1)
                    weights.append(weight)
                    valid_neighbors.append(neighbor)
            
            if not valid_neighbors:
                break
            
            # Weighted random selection
            weights = np.array(weights)
            probabilities = weights / weights.sum()
            
            next_node = np.random.choice(valid_neighbors, p=probabilities)
            used_edges.add((current_node, next_node))
            path.append(next_node)
            current_node = next_node
        
        return path if len(''.join(path)) >= target_length * 0.7 else None
    
    def _select_best_sequence(self, candidate_paths: List[List[str]], spectrum: List[str], 
                            target_length: int, k: int) -> str:
        """Select best sequence from candidates using multiple criteria"""
        best_sequence = None
        best_score = -1
        
        for path in candidate_paths:
            sequence = self._path_to_sequence(path, k)
            if not sequence:
                continue
            
            # Score based on multiple criteria
            spectrum_coverage = self._calculate_spectrum_coverage(sequence, spectrum, k)
            length_accuracy = 1.0 - abs(len(sequence) - target_length) / target_length
            nucleotide_balance = self._calculate_nucleotide_balance(sequence)
            
            combined_score = (spectrum_coverage * 0.5 + 
                            length_accuracy * 0.3 + 
                            nucleotide_balance * 0.2)
            
            if combined_score > best_score:
                best_score = combined_score
                best_sequence = sequence
        
        return best_sequence or self._enhanced_greedy_reconstruction(spectrum, target_length, k)
    
    def _calculate_spectrum_coverage(self, sequence: str, spectrum: List[str], k: int) -> float:
        """Calculate how well the sequence covers the original spectrum"""
        if len(sequence) < k:
            return 0.0
        
        sequence_kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
        sequence_kmer_set = set(sequence_kmers)
        spectrum_kmer_set = set(spectrum)
        
        if not spectrum_kmer_set:
            return 1.0
        
        intersection = len(sequence_kmer_set.intersection(spectrum_kmer_set))
        return intersection / len(spectrum_kmer_set)
    
    def _calculate_nucleotide_balance(self, sequence: str) -> float:
        """Calculate nucleotide composition balance (closer to 0.25 each is better)"""
        if not sequence:
            return 0.0
        
        counts = Counter(sequence)
        total = len(sequence)
        
        balance_score = 1.0
        for nucleotide in 'ACGT':
            freq = counts.get(nucleotide, 0) / total
            balance_score -= abs(freq - 0.25) * 0.5
        
        return max(0.0, balance_score)
    
    def _path_to_sequence(self, path: List[str], k: int) -> str:
        """Convert node path to DNA sequence"""
        if not path:
            return ""
        
        sequence = path[0]
        for i in range(1, len(path)):
            sequence += path[i][-1]
        
        return sequence
    
    def _enhanced_greedy_reconstruction(self, spectrum: List[str], target_length: int, k: int) -> str:
        """Enhanced greedy reconstruction with better heuristics"""
        if not spectrum:
            return 'A' * target_length
        
        # Start with most reliable k-mer
        kmer_scores = {}
        for kmer in set(spectrum):
            reliability = self._kmer_reliability.get(kmer, 0.5)
            kmer_scores[kmer] = reliability
        
        current_kmer = max(kmer_scores.keys(), key=lambda x: kmer_scores[x])
        result = current_kmer
        used = {current_kmer}
        
        # Build sequence with look-ahead
        while len(result) < target_length:
            suffix = result[-(k-1):]
            candidates = []
            
            for kmer in spectrum:
                if kmer not in used and kmer.startswith(suffix):
                    score = self._kmer_reliability.get(kmer, 0.5)
                    candidates.append((kmer, score))
            
            if not candidates:
                # Try with shorter suffix
                for overlap_len in range(k-2, 0, -1):
                    suffix = result[-overlap_len:]
                    for kmer in spectrum:
                        if kmer not in used and kmer.startswith(suffix):
                            score = self._kmer_reliability.get(kmer, 0.3)
                            candidates.append((kmer, score))
                    if candidates:
                        break
            
            if not candidates:
                # Random selection as last resort
                unused = [kmer for kmer in spectrum if kmer not in used]
                if unused:
                    next_kmer = max(unused, key=lambda x: self._kmer_reliability.get(x, 0.1))
                    result += next_kmer
                    used.add(next_kmer)
                else:
                    break
            else:
                # Select best candidate
                next_kmer = max(candidates, key=lambda x: x[1])[0]
                overlap_len = len(suffix)
                result += next_kmer[overlap_len:]
                used.add(next_kmer)
        
        return result[:target_length] 