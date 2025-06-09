import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from typing import List, Dict
import os

class BenchmarkVisualizer:
    def __init__(self, output_dir: str = "benchmark_results/plots"):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # Set style
        plt.style.use('seaborn-v0_8')  # Using a specific seaborn style version
        sns.set_theme()  # Set seaborn theme
    
    def plot_accuracy_by_k(self, results: List[Dict]):
        """Plot accuracy vs k-mer size."""
        df = pd.DataFrame(results)
        
        plt.figure(figsize=(10, 6))
        sns.boxplot(x='k', y='coverage', data=df)
        plt.title('Accuracy by K-mer Size')
        plt.xlabel('K-mer Size')
        plt.ylabel('Spectrum Coverage (%)')
        plt.grid(True, alpha=0.3)
        
        # Save plot
        plt.savefig(os.path.join(self.output_dir, 'accuracy_by_k.png'))
        plt.close()
    
    def plot_accuracy_by_length(self, results: List[Dict]):
        """Plot accuracy vs sequence length."""
        df = pd.DataFrame(results)
        
        plt.figure(figsize=(10, 6))
        sns.scatterplot(x='sequence_length', y='coverage', data=df, hue='k')
        plt.title('Accuracy by Sequence Length')
        plt.xlabel('Sequence Length')
        plt.ylabel('Spectrum Coverage (%)')
        plt.grid(True, alpha=0.3)
        
        # Save plot
        plt.savefig(os.path.join(self.output_dir, 'accuracy_by_length.png'))
        plt.close()
    
    def plot_runtime_by_k(self, results: List[Dict]):
        """Plot runtime vs k-mer size."""
        df = pd.DataFrame(results)
        
        plt.figure(figsize=(10, 6))
        sns.boxplot(x='k', y='runtime', data=df)
        plt.title('Runtime by K-mer Size')
        plt.xlabel('K-mer Size')
        plt.ylabel('Runtime (seconds)')
        plt.grid(True, alpha=0.3)
        
        # Save plot
        plt.savefig(os.path.join(self.output_dir, 'runtime_by_k.png'))
        plt.close()
    
    def plot_runtime_by_length(self, results: List[Dict]):
        """Plot runtime vs sequence length."""
        df = pd.DataFrame(results)
        
        plt.figure(figsize=(10, 6))
        sns.scatterplot(x='sequence_length', y='runtime', data=df, hue='k')
        plt.title('Runtime by Sequence Length')
        plt.xlabel('Sequence Length')
        plt.ylabel('Runtime (seconds)')
        plt.grid(True, alpha=0.3)
        
        # Save plot
        plt.savefig(os.path.join(self.output_dir, 'runtime_by_length.png'))
        plt.close()
    
    def plot_all_metrics(self, results: List[Dict]):
        """Generate all plots."""
        self.plot_accuracy_by_k(results)
        self.plot_accuracy_by_length(results)
        self.plot_runtime_by_k(results)
        self.plot_runtime_by_length(results)
        
        # Create summary statistics
        df = pd.DataFrame(results)
        summary = df.groupby('k').agg({
            'coverage': ['mean', 'std', 'min', 'max'],
            'runtime': ['mean', 'std', 'min', 'max']
        }).round(2)
        
        # Save summary to CSV
        summary.to_csv(os.path.join(self.output_dir, 'summary_statistics.csv')) 