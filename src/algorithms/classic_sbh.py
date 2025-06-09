from typing import List, Dict, Set, Tuple, Optional
import networkx as nx
from collections import defaultdict, Counter

class ClassicSBH:
    """Classic Sequencing by Hybridization algorithm implementation"""
    
    def __init__(self):
        self._graph = None
    
    def reconstruct(self, spectrum: List[str], target_length: int, k: int) -> str:
        """
        Reconstruct DNA sequence from its spectrum using classic SBH algorithm
        
        Args:
            spectrum: List of k-mers from hybridization
            target_length: Expected length of the DNA sequence
            k: Length of oligonucleotides in spectrum
            
        Returns:
            Reconstructed DNA sequence
            
        Raises:
            ValueError: If input parameters are invalid
        """
        # Validate input
        if not spectrum:
            raise ValueError("Spectrum cannot be empty")
        if any(len(kmer) != k for kmer in spectrum):
            raise ValueError("All k-mers must have length k")
        if target_length < k:
            raise ValueError("Target length must be at least k")
        if k < 1:
            raise ValueError("k must be positive")
            
        print("\n=== Starting DNA Sequence Reconstruction ===")
        print(f"Target length: {target_length}")
        print(f"K-mer length: {k}")
        print(f"Spectrum size: {len(spectrum)}")
        print(f"Unique k-mers: {len(set(spectrum))}")
        
        try:
            print("\nDebug: Building de Bruijn graph...")
            self._graph = self._build_graph(spectrum, k)
            print(f"Debug: Graph built with {self._graph.number_of_nodes()} nodes and {self._graph.number_of_edges()} edges")
            
            try:
                print("\nDebug: Attempting to find Eulerian path...")
                path = self._find_eulerian_path(self._graph)
                if path:
                    sequence = self._sequence_from_path(path, k)
                    if len(sequence) != target_length:
                        print(f"Warning: Sequence length mismatch. Expected {target_length}, got {len(sequence)}")
                    return sequence
                else:
                    print("Debug: No valid Eulerian path found")
                    return self._greedy_reconstruction(spectrum, target_length, k)
                
            except Exception as e:
                print(f"\nError during reconstruction: {e}")
                print("Debug: Attempting greedy reconstruction as fallback")
                return self._greedy_reconstruction(spectrum, target_length, k)
                
        except Exception as e:
            print(f"\nError during reconstruction: {e}")
            print("Debug: Attempting greedy reconstruction as fallback")
            return self._greedy_reconstruction(spectrum, target_length, k)
    
    def _build_graph(self, spectrum: List[str], k: int) -> nx.MultiDiGraph:
        """Build de Bruijn graph from spectrum"""
        # Count k-mer occurrences
        kmer_counts = Counter(spectrum)
        
        # Create a multigraph to handle duplicate k-mers
        graph = nx.MultiDiGraph()
        
        # Add edges for each k-mer occurrence
        for kmer, count in kmer_counts.items():
            prefix = kmer[:-1]
            suffix = kmer[1:]
            # Add edges with multiplicity
            for _ in range(count):
                graph.add_edge(prefix, suffix, kmer=kmer)
            
        return graph
    
    def _find_eulerian_path(self, graph: nx.DiGraph) -> Optional[List[Tuple[str, str, int]]]:
        """Find Eulerian path in the de Bruijn graph"""
        if not graph.edges:
            return None
            
        # Find connected components
        components = list(nx.weakly_connected_components(graph))
        print(f"Debug: Found {len(components)} connected components")
        
        if len(components) > 1:
            # Try to find paths in each component and connect them
            component_paths = []
            for comp in components:
                subgraph = graph.subgraph(comp).copy()
                path = self._find_path_in_component(subgraph)
                if path:
                    component_paths.append(path)
                    
            if not component_paths:
                print("Debug: No valid paths found in components")
                return None
                
            # Try to connect paths using overlaps
            final_path = []
            used_paths = set()
            current_path = component_paths[0]
            used_paths.add(0)
            final_path.extend(current_path)
            
            while len(used_paths) < len(component_paths):
                best_overlap = -1
                best_path_idx = -1
                best_path = None
                
                last_kmer = self._get_kmer_from_path(final_path[-1])
                
                for i, path in enumerate(component_paths):
                    if i in used_paths:
                        continue
                        
                    first_kmer = self._get_kmer_from_path(path[0])
                    overlap = self._find_overlap(last_kmer, first_kmer)
                    
                    if overlap > best_overlap:
                        best_overlap = overlap
                        best_path_idx = i
                        best_path = path
                        
                if best_path is None:
                    print("Debug: Could not connect all components")
                    break
                    
                used_paths.add(best_path_idx)
                final_path.extend(best_path)
                
            return final_path
            
        # Single component - try to find Eulerian path
        return self._find_path_in_component(graph)
        
    def _find_path_in_component(self, graph: nx.DiGraph) -> Optional[List[Tuple[str, str, int]]]:
        """Find path in a single component"""
        # Find start node (node with in_degree < out_degree)
        start_node = None
        for node in graph.nodes():
            in_deg = graph.in_degree(node)
            out_deg = graph.out_degree(node)
            if out_deg > in_deg:
                start_node = node
                break
                
        # If no start node found, try any node with outgoing edges
        if not start_node:
            for node in graph.nodes():
                if graph.out_degree(node) > 0:
                    start_node = node
                    break
                    
        if not start_node:
            print("Debug: No valid start node found")
            return None
            
        print(f"Debug: Found start node: {start_node}")
        
        # Try to find path using DFS with backtracking
        path = []
        used_edges = set()
        
        def dfs(node: str) -> bool:
            if len(used_edges) == graph.number_of_edges():
                return True
                
            # Sort edges by their weight (if available) or degree
            edges = list(graph.out_edges(node, keys=True))
            edges.sort(key=lambda e: (
                -graph[e[0]][e[1]][e[2]].get('weight', 0),
                -graph.out_degree(e[1])
            ))
            
            for edge in edges:
                if edge not in used_edges:
                    used_edges.add(edge)
                    path.append(edge)
                    
                    if dfs(edge[1]):
                        return True
                        
                    used_edges.remove(edge)
                    path.pop()
                    
            return False
            
        if dfs(start_node):
            print(f"Debug: Found valid path of length {len(path)}")
            return path
            
        print("Debug: No valid path found")
        return None
        
    def _get_kmer_from_path(self, edge: Tuple[str, str, int]) -> str:
        """Get k-mer from path edge"""
        return self._graph[edge[0]][edge[1]][edge[2]]['kmer']
        
    def _find_overlap(self, kmer1: str, kmer2: str) -> int:
        """Find overlap length between two k-mers"""
        k = len(kmer1)
        for i in range(k-1, 0, -1):
            if kmer1.endswith(kmer2[:i]):
                return i
        return 0
    
    def _sequence_from_path(self, path: List[Tuple[str, str, int]], k: int) -> str:
        """Reconstruct DNA sequence from Eulerian path"""
        if not path:
            return ""
            
        print("\nDebug: Reconstructing sequence from path...")
        # Get all k-mers in the path
        kmers = [self._get_kmer_from_edge(edge) for edge in path]
        print(f"Debug: K-mers in path: {kmers}")
        
        # Start with the first k-mer
        sequence = list(kmers[0])
        print(f"Debug: Starting sequence: {''.join(sequence)}")
        
        # Add one nucleotide at a time from each subsequent k-mer
        for i in range(1, len(kmers)):
            # Verify overlap
            if kmers[i][:-1] != kmers[i-1][1:]:
                raise ValueError(f"K-mers {kmers[i-1]} and {kmers[i]} do not overlap properly")
            sequence.append(kmers[i][-1])
            print(f"Debug: After adding {kmers[i][-1]}: {''.join(sequence)}")
            
        return ''.join(sequence)
    
    def _get_kmer_from_edge(self, edge: Tuple[str, str, int]) -> str:
        """Get the k-mer associated with an edge"""
        u, v, key = edge  # Get source, target nodes and edge key
        return self._graph[u][v][key]['kmer']
    
    def _find_best_kmer_match(self, current_suffix: str, spectrum: List[str], used: Set[str], k: int) -> Tuple[str, float]:
        """Find the best matching k-mer from the spectrum"""
        best_match = None
        best_score = -1
        
        # Convert suffix to string if it's a list
        if isinstance(current_suffix, list):
            current_suffix = ''.join(current_suffix)
            
        # First try exact suffix matches
        suffix_len = min(len(current_suffix), k-1)
        suffix = current_suffix[-suffix_len:]
        
        # Look for exact matches first
        exact_matches = []
        for kmer in spectrum:
            if kmer in used:
                continue
            if kmer.startswith(suffix):
                # Calculate score based on frequency and overlap length
                freq = spectrum.count(kmer)
                score = len(suffix) + 0.2 * freq
                exact_matches.append((kmer, score))
                
        if exact_matches:
            # Sort by score and frequency
            exact_matches.sort(key=lambda x: (-x[1], -spectrum.count(x[0])))
            return exact_matches[0]
            
        # If no exact matches, try partial overlaps with higher weight on longer overlaps
        for kmer in spectrum:
            if kmer in used:
                continue
                
            # Find longest overlap
            for i in range(min(len(current_suffix), k-1), 2, -1):  # Require at least 3 bp overlap
                if current_suffix[-i:] == kmer[:i]:
                    # Calculate score based on overlap length, frequency, and position
                    freq = spectrum.count(kmer)
                    position_score = 1.0 if i == k-1 else 0.5  # Higher score for full k-1 overlap
                    score = i + 0.2 * freq + position_score
                    
                    if score > best_score:
                        best_score = score
                        best_match = kmer
                        break
                        
        return best_match, best_score

    def _greedy_reconstruction(self, spectrum: List[str], target_length: int, k: int) -> str:
        """Greedy reconstruction when no Eulerian path exists"""
        print("\n=== Starting Greedy Reconstruction ===")
        
        import time
        start_time = time.time()
        max_runtime = 30  # Maximum 30 seconds for reconstruction
        
        # Calculate k-mer frequencies and overlap scores
        kmer_scores = {}
        for kmer in spectrum:
            score = 0
            for other in spectrum:
                if kmer != other:
                    # Check both prefix and suffix overlaps
                    for i in range(k-1, 2, -1):
                        if kmer.endswith(other[:i]) or kmer.startswith(other[-i:]):
                            score += i * spectrum.count(other)  # Weight by frequency
            kmer_scores[kmer] = score + 0.2 * spectrum.count(kmer)
            
        # Start with most promising k-mer
        current = max(spectrum, key=lambda x: kmer_scores[x])
        result = current
        used = {current}
        
        print(f"Debug: Starting with best overlapping k-mer: {current} (score: {kmer_scores[current]:.1f})")
        print(f"\nDebug: Current sequence length: {len(result)}/{target_length}")
        
        consecutive_failures = 0
        max_consecutive_failures = 5
        backtrack_points = []  # Store points where we can backtrack
        failed_paths = set()  # Track failed path states to prevent infinite loops
        min_overlap = 3  # Minimum required overlap
        max_backtrack_attempts = 10  # Limit total backtrack attempts
        backtrack_count = 0
        iteration_count = 0
        max_iterations = target_length * 5  # Prevent infinite loops
        
        while len(result) < target_length and backtrack_count < max_backtrack_attempts:
            iteration_count += 1
            
            # Check timeout and iteration limits
            if time.time() - start_time > max_runtime:
                print(f"Debug: Reconstruction timeout after {max_runtime} seconds")
                break
            if iteration_count > max_iterations:
                print(f"Debug: Reached maximum iterations ({max_iterations}), stopping")
                break
            # Get suffix for matching
            suffix = result[-k+1:] if len(result) >= k-1 else result
            
            # Create a state key to detect repeated failures
            state_key = (suffix, frozenset(used), len(result))
            
            # Check if we've already failed from this state
            if state_key in failed_paths:
                print(f"Debug: Detected previously failed state, forcing backtrack")
                consecutive_failures = max_consecutive_failures
            else:
                print(f"Debug: Looking for k-mer with prefix {suffix}")
                
                # Find best matching k-mer
                next_kmer, score = self._find_best_kmer_match(suffix, spectrum, used, k)
                
                if next_kmer and score > min_overlap:
                    print(f"Debug: Found match {next_kmer} with score {score:.1f}")
                    used.add(next_kmer)
                    result = result + next_kmer[-1]  # Add only the last nucleotide
                    print(f"Debug: Current sequence: {result} (length: {len(result)}/{target_length})")
                    consecutive_failures = 0
                    
                    # Store good points for backtracking (limit the number stored)
                    if score >= k/2 and len(backtrack_points) < 20:
                        backtrack_points.append((len(result), set(used), result))
                    continue
                else:
                    consecutive_failures += 1
            
            # Handle failures and backtracking
            if consecutive_failures >= max_consecutive_failures:
                print(f"Debug: Failed to find good matches after {consecutive_failures} attempts")
                
                # Mark this state as failed
                failed_paths.add(state_key)
                
                # Try to backtrack
                if backtrack_points:
                    backtrack_count += 1
                    pos, used_kmers, seq = backtrack_points.pop()
                    print(f"Debug: Backtracking to position {pos} (attempt {backtrack_count})")
                    result = seq
                    used = used_kmers.copy()  # Make a copy to avoid reference issues
                    consecutive_failures = 0
                    
                    # Remove some recent failed paths to allow exploration
                    if backtrack_count > 5:
                        failed_paths.clear()
                        print("Debug: Cleared failed paths cache to enable new exploration")
                    continue
                else:
                    print("Debug: No backtrack points available")
                    break
                    
            # Try to recover by using most promising unused k-mer
            unused_kmers = [(k, kmer_scores[k]) for k in spectrum if k not in used]
            if unused_kmers:
                next_kmer = max(unused_kmers, key=lambda x: x[1])[0]
                print(f"Debug: Recovering with high-scoring unused k-mer {next_kmer}")
                used.add(next_kmer)
                result = result + next_kmer[-1]
                consecutive_failures = 0
            else:
                print("Debug: No unused k-mers left")
                break
                
        if backtrack_count >= max_backtrack_attempts:
            print(f"Debug: Reached maximum backtrack attempts ({max_backtrack_attempts}), stopping")
                    
        return result 