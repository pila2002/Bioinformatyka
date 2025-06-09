#!/usr/bin/env python3
"""
Script to generate professional charts for the SBH project report
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path

# Set style for professional look
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

def create_charts_directory():
    """Create directory for charts if it doesn't exist"""
    charts_dir = Path("report_charts")
    charts_dir.mkdir(exist_ok=True)
    return charts_dir

def generate_candidate_size_chart(charts_dir):
    """Generate candidate_size performance chart"""
    # Data from our tests
    candidate_sizes = [3, 5, 8, 10, 15, 20, 25, 30]
    avg_accuracy = [49.11, 47.11, 51.11, 47.00, 50.78, 51.22, 50.44, 48.78]
    avg_time = [0.107, 0.025, 0.043, 0.038, 0.020, 0.014, 0.086, 0.024]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Accuracy chart
    bars1 = ax1.bar(candidate_sizes, avg_accuracy, color='steelblue', alpha=0.7)
    ax1.set_xlabel('Rozmiar candidate_size', fontsize=12)
    ax1.set_ylabel('Średnia dokładność (%)', fontsize=12)
    ax1.set_title('Wpływ parametru candidate_size na dokładność', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Highlight optimal value
    optimal_idx = candidate_sizes.index(20)
    bars1[optimal_idx].set_color('red')
    bars1[optimal_idx].set_alpha(0.8)
    
    # Add value labels
    for i, v in enumerate(avg_accuracy):
        ax1.text(candidate_sizes[i], v + 0.2, f'{v:.1f}%', ha='center', fontweight='bold')
    
    # Time chart
    bars2 = ax2.bar(candidate_sizes, avg_time, color='darkgreen', alpha=0.7)
    ax2.set_xlabel('Rozmiar candidate_size', fontsize=12)
    ax2.set_ylabel('Średni czas wykonania (s)', fontsize=12)
    ax2.set_title('Wpływ parametru candidate_size na czas wykonania', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # Highlight optimal value
    bars2[optimal_idx].set_color('red')
    bars2[optimal_idx].set_alpha(0.8)
    
    # Add value labels
    for i, v in enumerate(avg_time):
        ax2.text(candidate_sizes[i], v + max(avg_time)*0.02, f'{v:.3f}s', ha='center', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(charts_dir / 'candidate_size_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✅ Wykres candidate_size zapisany")

def generate_dna_length_chart(charts_dir):
    """Generate DNA length impact chart"""
    # Data from our tests
    lengths = [200, 400, 600]
    accuracies = [46.10, 48.60, 50.23]
    times = [0.056, 0.068, 0.097]
    variances = [3.28, 2.77, 0.64]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Accuracy and variance
    ax1_twin = ax1.twinx()
    
    line1 = ax1.plot(lengths, accuracies, 'o-', color='blue', linewidth=3, markersize=8, label='Dokładność')
    line2 = ax1_twin.plot(lengths, variances, 's--', color='red', linewidth=2, markersize=6, label='Wariancja')
    
    ax1.set_xlabel('Długość sekwencji DNA (nukleotydów)', fontsize=12)
    ax1.set_ylabel('Średnia dokładność (%)', fontsize=12, color='blue')
    ax1_twin.set_ylabel('Wariancja (%)', fontsize=12, color='red')
    ax1.set_title('Wpływ długości sekwencji na jakość rekonstrukcji', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Add value labels
    for i, (x, y) in enumerate(zip(lengths, accuracies)):
        ax1.text(x, y + 0.3, f'{y:.1f}%', ha='center', fontweight='bold', color='blue')
    
    for i, (x, y) in enumerate(zip(lengths, variances)):
        ax1_twin.text(x, y + 0.1, f'{y:.1f}%', ha='center', fontweight='bold', color='red')
    
    # Combine legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1_twin.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='center right')
    
    # Time chart
    bars = ax2.bar(lengths, times, color='orange', alpha=0.7, width=50)
    ax2.set_xlabel('Długość sekwencji DNA (nukleotydów)', fontsize=12)
    ax2.set_ylabel('Średni czas wykonania (s)', fontsize=12)
    ax2.set_title('Wpływ długości sekwencji na czas wykonania', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # Add value labels
    for i, v in enumerate(times):
        ax2.text(lengths[i], v + max(times)*0.02, f'{v:.3f}s', ha='center', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(charts_dir / 'dna_length_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✅ Wykres długości DNA zapisany")

def generate_kmer_size_chart(charts_dir):
    """Generate k-mer size impact chart"""
    # Data from our tests
    k_sizes = [7, 8, 9]
    accuracies = [46.75, 48.60, 44.20]
    times = [0.356, 0.068, 1.038]
    contigs = [38, 35, 0]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Accuracy and contigs
    ax1_twin = ax1.twinx()
    
    bars1 = ax1.bar([x-0.2 for x in k_sizes], accuracies, width=0.4, color='purple', alpha=0.7, label='Dokładność')
    bars2 = ax1_twin.bar([x+0.2 for x in k_sizes], contigs, width=0.4, color='green', alpha=0.7, label='Liczba kontigów')
    
    ax1.set_xlabel('Rozmiar k-meru', fontsize=12)
    ax1.set_ylabel('Średnia dokładność (%)', fontsize=12, color='purple')
    ax1_twin.set_ylabel('Liczba kontigów w fazie 1', fontsize=12, color='green')
    ax1.set_title('Wpływ rozmiaru k-meru na jakość rekonstrukcji', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(k_sizes)
    
    # Add value labels
    for i, v in enumerate(accuracies):
        ax1.text(k_sizes[i]-0.2, v + 0.5, f'{v:.1f}%', ha='center', fontweight='bold', color='purple')
    
    for i, v in enumerate(contigs):
        ax1_twin.text(k_sizes[i]+0.2, v + 1, f'{v}', ha='center', fontweight='bold', color='green')
    
    # Combine legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1_twin.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
    
    # Time chart (log scale due to large differences)
    bars = ax2.bar(k_sizes, times, color='brown', alpha=0.7)
    ax2.set_xlabel('Rozmiar k-meru', fontsize=12)
    ax2.set_ylabel('Średni czas wykonania (s)', fontsize=12)
    ax2.set_title('Wpływ rozmiaru k-meru na czas wykonania', fontsize=14, fontweight='bold')
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3)
    ax2.set_xticks(k_sizes)
    
    # Add value labels
    for i, v in enumerate(times):
        ax2.text(k_sizes[i], v * 1.2, f'{v:.3f}s', ha='center', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(charts_dir / 'kmer_size_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✅ Wykres rozmiaru k-meru zapisany")

def generate_error_impact_chart(charts_dir):
    """Generate error level impact chart"""
    # Data from our tests
    error_levels = ['2%+2%', '5%+5%', '10%+10%']
    accuracies = [49.00, 48.60, 49.35]
    times = [0.109, 0.068, 0.039]
    stability = [50.00-46.25, 52.50-45.50, 50.25-48.50]  # Range as stability measure
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Accuracy with error bars (stability)
    bars1 = ax1.bar(error_levels, accuracies, color='teal', alpha=0.7)
    ax1.set_xlabel('Poziom błędów (pozytywne + negatywne)', fontsize=12)
    ax1.set_ylabel('Średnia dokładność (%)', fontsize=12)
    ax1.set_title('Wpływ poziomu błędów na dokładność rekonstrukcji', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Add error bars representing stability (smaller = more stable)
    ax1.errorbar(range(len(error_levels)), accuracies, yerr=[[s/2 for s in stability], [s/2 for s in stability]], 
                 fmt='none', color='red', capsize=5, capthick=2)
    
    # Add value labels
    for i, v in enumerate(accuracies):
        ax1.text(i, v + 0.8, f'{v:.1f}%', ha='center', fontweight='bold')
    
    # Time chart
    bars2 = ax2.bar(error_levels, times, color='coral', alpha=0.7)
    ax2.set_xlabel('Poziom błędów (pozytywne + negatywne)', fontsize=12)
    ax2.set_ylabel('Średni czas wykonania (s)', fontsize=12)
    ax2.set_title('Wpływ poziomu błędów na czas wykonania', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # Add value labels
    for i, v in enumerate(times):
        ax2.text(i, v + max(times)*0.02, f'{v:.3f}s', ha='center', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(charts_dir / 'error_impact_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✅ Wykres wpływu błędów zapisany")

def generate_algorithm_comparison_chart(charts_dir):
    """Generate algorithm comparison chart"""
    algorithms = ['Klasyczny SBH', 'Trzyfazowy SBH']
    accuracies = [24.1, 24.2]
    times = [0.400, 0.101]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    # Accuracy comparison
    bars1 = ax1.bar(algorithms, accuracies, color=['lightcoral', 'lightgreen'], alpha=0.8)
    ax1.set_ylabel('Średnia dokładność (%)', fontsize=12)
    ax1.set_title('Porównanie dokładności algorytmów', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Add value labels
    for i, v in enumerate(accuracies):
        ax1.text(i, v + 0.2, f'{v:.1f}%', ha='center', fontweight='bold')
    
    # Time comparison
    bars2 = ax2.bar(algorithms, times, color=['lightcoral', 'lightgreen'], alpha=0.8)
    ax2.set_ylabel('Średni czas wykonania (s)', fontsize=12)
    ax2.set_title('Porównanie wydajności algorytmów', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # Add value labels and improvement annotation
    for i, v in enumerate(times):
        ax2.text(i, v + max(times)*0.02, f'{v:.3f}s', ha='center', fontweight='bold')
    
    # Add improvement annotation
    improvement = times[0] / times[1]
    ax2.annotate(f'{improvement:.1f}x szybszy!', 
                xy=(1, times[1]), xytext=(0.5, times[0]*0.7),
                arrowprops=dict(arrowstyle='->', color='red', lw=2),
                fontsize=12, fontweight='bold', color='red',
                ha='center')
    
    plt.tight_layout()
    plt.savefig(charts_dir / 'algorithm_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✅ Wykres porównania algorytmów zapisany")

def main():
    """Generate all charts for the report"""
    print("🎨 Generowanie wykresów do sprawozdania...")
    
    charts_dir = create_charts_directory()
    
    generate_candidate_size_chart(charts_dir)
    generate_dna_length_chart(charts_dir)
    generate_kmer_size_chart(charts_dir)
    generate_error_impact_chart(charts_dir)
    generate_algorithm_comparison_chart(charts_dir)
    
    print(f"\n✅ Wszystkie wykresy zostały zapisane w katalogu: {charts_dir}")
    print("📊 Gotowe wykresy:")
    for chart_file in sorted(charts_dir.glob("*.png")):
        print(f"   - {chart_file.name}")

if __name__ == "__main__":
    main() 