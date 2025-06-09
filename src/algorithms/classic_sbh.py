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
                
            for _, next_node, key in graph.out_edges(node, keys=True):
                edge = (node, next_node, key)
                if edge not in used_edges:
                    used_edges.add(edge)
                    path.append(edge)
                    
                    if dfs(next_node):
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
        # Create a copy of the graph to mark edges
        working_graph = graph.copy()
        for u, v, k in working_graph.edges(keys=True):
            working_graph[u][v][k]['used'] = False
        
        # Try to find Eulerian path
        try:
            path = []
            current = start_node
            stack = [(current, [])]  # (node, path_so_far)
            
            while stack:
                current, current_path = stack.pop()
                
                # Get unused outgoing edges
                out_edges = [(u, v, k) for u, v, k in working_graph.out_edges(current, keys=True)
                           if not working_graph[u][v][k]['used']]
                
                if not out_edges:
                    if len(current_path) > len(path):
                        path = current_path
                    continue
                    
                # Try each unused edge
                for edge in out_edges:
                    working_graph[edge[0]][edge[1]][edge[2]]['used'] = True
                    new_path = current_path + [edge]
                    stack.append((edge[1], new_path))
            
            # Validate path
            if not path:
                print("Debug: No path found")
                return None
                
            # Check if path uses all edges
            edge_count = graph.number_of_edges()
            if len(path) != edge_count:
                print(f"Debug: Path does not use all edges: found {len(path)} of {edge_count}")
                return None
                
            # Check if path is connected
            nodes_in_path = {start_node} | {edge[1] for edge in path}
            if len(nodes_in_path) != graph.number_of_nodes():
                print(f"Debug: Path does not visit all nodes: visited {len(nodes_in_path)} of {graph.number_of_nodes()}")
                return None
                
            print(f"Debug: Found valid Eulerian path of length {len(path)}")
            return path
            
        except Exception as e:
            print(f"Debug: Error finding Eulerian path: {str(e)}")
            return None
    
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
                score = len(suffix) + 0.1 * spectrum.count(kmer)
                exact_matches.append((kmer, score))
                
        if exact_matches:
            # Sort by score and frequency
            exact_matches.sort(key=lambda x: (-x[1], -spectrum.count(x[0])))
            return exact_matches[0]
            
        # If no exact matches, try partial overlaps
        for kmer in spectrum:
            if kmer in used:
                continue
                
            # Find longest overlap
            for i in range(min(len(current_suffix), k-1), 2, -1):  # Require at least 3 bp overlap
                if current_suffix[-i:] == kmer[:i]:
                    score = i + 0.1 * spectrum.count(kmer)
                    if score > best_score:
                        best_score = score
                        best_match = kmer
                        break
                        
        return best_match, best_score

    def _greedy_reconstruction(self, spectrum: List[str], target_length: int, k: int) -> str:
        """Greedy reconstruction when no Eulerian path exists"""
        print("\n=== Starting Greedy Reconstruction ===")
        
        # Start with most frequent k-mer that has most overlaps
        kmer_scores = {}
        for kmer in spectrum:
            score = 0
            for other in spectrum:
                if kmer != other:
                    # Check both prefix and suffix overlaps
                    for i in range(k-1, 2, -1):
                        if kmer.endswith(other[:i]) or kmer.startswith(other[-i:]):
                            score += i
                            break
            kmer_scores[kmer] = score + 0.1 * spectrum.count(kmer)
            
        current = max(spectrum, key=lambda x: kmer_scores[x])
        result = current
        used = {current}
        
        print(f"Debug: Starting with best overlapping k-mer: {current} (score: {kmer_scores[current]:.1f})")
        print(f"\nDebug: Current sequence length: {len(result)}/{target_length}")
        
        consecutive_failures = 0
        max_consecutive_failures = 5
        backtrack_points = []  # Store points where we can backtrack
        
        while len(result) < target_length:
            # Get suffix for matching
            suffix = result[-k+1:] if len(result) >= k-1 else result
            print(f"Debug: Looking for k-mer with prefix {suffix}")
            
            # Find best matching k-mer
            next_kmer, score = self._find_best_kmer_match(suffix, spectrum, used, k)
            
            if next_kmer and score > 3:  # Require at least 3 bp overlap
                print(f"Debug: Found match {next_kmer} with score {score:.1f}")
                used.add(next_kmer)
                result = result + next_kmer[-1]  # Add only the last nucleotide
                print(f"Debug: Current sequence: {result} (length: {len(result)}/{target_length})")
                consecutive_failures = 0
                
                # Store good points for backtracking
                if score >= k/2:
                    backtrack_points.append((len(result), set(used), result))
            else:
                consecutive_failures += 1
                if consecutive_failures >= max_consecutive_failures:
                    print(f"Debug: Failed to find good matches after {consecutive_failures} attempts")
                    
                    # Try to backtrack
                    if backtrack_points:
                        pos, used_kmers, seq = backtrack_points.pop()
                        print(f"Debug: Backtracking to position {pos}")
                        result = seq
                        used = used_kmers
                        consecutive_failures = 0
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
                else:
                    print("Debug: No unused k-mers left")
                    break
                    
        return result 