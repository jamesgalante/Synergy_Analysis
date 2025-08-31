import hicstraw
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os

def visualize_hic(hic_file, chrom, start, end, resolution, 
                  save_path=None, title=None, log_transform=True, markers=None):
    # Clear any existing plots to avoid multiple colorbars
    plt.clf()
    
    # Open Hi-C file
    hic = hicstraw.HiCFile(hic_file)
    
    # Get matrix data
    matrix = hic.getMatrixZoomData(chrom, chrom, "observed", "SCALE", "BP", resolution)
    
    # Get contact records for the region
    records = matrix.getRecords(start, end, start, end)
    
    # Build contact matrix
    size = (end - start) // resolution
    contact_matrix = np.zeros((size, size))
    
    # Fill in the contact matrix
    for record in records:
        i = (record.binX - start) // resolution
        j = (record.binY - start) // resolution
        if 0 <= i < size and 0 <= j < size:
            contact_matrix[i, j] = record.counts
            contact_matrix[j, i] = record.counts  # Make symmetric
    
    # Apply log transformation if requested
    if log_transform:
        plot_matrix = np.log1p(contact_matrix)
        cbar_label = 'log(contacts + 1)'
    else:
        plot_matrix = contact_matrix
        cbar_label = 'contacts'
    
    # Create the plot
    plt.figure(figsize=(12, 8))
    
    # Create heatmap (SINGLE call - this was duplicated before!)
    ax = sns.heatmap(plot_matrix, 
                     cmap='Reds', 
                     square=True,
                     cbar_kws={'label': cbar_label, 'shrink': 0.8},
                     xticklabels=False, 
                     yticklabels=False)
    
    # Add markers if provided
    if markers is not None:
        # Convert genomic coordinates to matrix indices
        def pos_to_matrix_idx(genomic_pos):
            if start <= genomic_pos < end:
                return (genomic_pos - start) // resolution
            return None
        
        # Process each marker
        for label, marker_info in markers.items():
            # Handle different input formats
            if isinstance(marker_info, (int, float)):
                # Simple format: just position
                pos = marker_info
                color = 'blue'
                marker_style = 'o'
            elif isinstance(marker_info, tuple):
                # Tuple format: (pos, color, marker)
                pos, color, marker_style = marker_info
            elif isinstance(marker_info, dict):
                # Dict format: {'pos': x, 'color': 'blue', 'marker': 'o'}
                pos = marker_info['pos']
                color = marker_info.get('color', 'blue')
                marker_style = marker_info.get('marker', 'o')
            else:
                continue
            
            # Convert to matrix coordinates
            idx = pos_to_matrix_idx(pos)
            if idx is not None:
                # Add crosshairs
                ax.axhline(y=idx + 0.5, color=color, linestyle='-', linewidth=2, alpha=0.8, label=label)
                ax.axvline(x=idx + 0.5, color=color, linestyle='-', linewidth=2, alpha=0.8)
                # Add point marker
                ax.scatter(idx + 0.5, idx + 0.5, color=color, s=100, marker=marker_style, 
                          edgecolor='white', linewidth=1)
        
        # Add legend positioned to avoid colorbar
        ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1), frameon=True, fancybox=True, shadow=True)
    
    # Set title
    if title is None:
        title = f'Hi-C Contacts: {chrom}:{start:,}-{end:,} (res={resolution:,}bp)'
    plt.title(title, fontsize=14)
    
    # Add position labels
    plt.xlabel(f'Position ({chrom})', fontsize=12)
    plt.ylabel(f'Position ({chrom})', fontsize=12)
    
    # Save or show
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Plot saved to: {save_path}")
    else:
        plt.show()
    
    return contact_matrix

# Read synergy predictions
synergy_df = pd.read_csv(snakemake.input.synergy_predictions, sep='\t')

# Get parameters from Snakemake
hic_file = snakemake.input.hic
hic_res = snakemake.params.hic_res
hic_res_small = snakemake.params.hic_res_small  # Default to hic_res if not provided
plot_window = snakemake.params.plot_window

# Debug: Print the actual parameter values
print(f">>> DEBUG: hic_res = {hic_res} (type: {type(hic_res)})")
print(f">>> DEBUG: hic_res_small = {hic_res_small} (type: {type(hic_res_small)})")
print(f">>> DEBUG: plot_window = {plot_window} (type: {type(plot_window)})")

# Create output directory for plots
plot_dir = snakemake.params.plots_dir
os.makedirs(plot_dir, exist_ok=True)

print(f"Processing {len(synergy_df)} synergistic pairs...")
print(f"Hi-C resolution: {hic_res}, small resolution: {hic_res_small}")
print(f"Plot window: {plot_window}")

# Create plots for each synergistic pair
plot_count = 0
for idx, row in synergy_df.iterrows():
    try:
        # Get chromosomes (assuming both elements are on same chromosome)
        chrom = row['#chr_element1']
        
        # Calculate average of element1 start and TSS location
        element1_start = row['start_element1']
        tss_pos = row['TargetGeneTSS_element1']
        region_center = (element1_start + tss_pos) // 2
        
        # Calculate distance between element1 start and TSS
        distance_e1_tss = abs(element1_start - tss_pos)
        
        # Calculate minimum window size (100 resolution bins) using ORIGINAL resolution
        min_window = 100 * hic_res
        
        # Debug: Print calculations for first few pairs
        if idx < 3:
            print(f">>> DEBUG pair {idx}: hic_res={hic_res}, min_window={min_window}")
        
        # Error check: minimum window can't exceed plotting window
        if min_window > plot_window:
            raise ValueError(f"Minimum window size ({min_window:,} bp) exceeds maximum plotting window ({plot_window:,} bp). Reduce resolution or increase plotting window.")
        
        # Calculate adaptive window: 2x distance on each side of average = 4x total distance
        adaptive_window = 4 * distance_e1_tss
        
        # Determine if we need to use smaller resolution
        use_small_resolution = adaptive_window < min_window
        
        # Apply window size constraints
        window_size = max(adaptive_window, min_window)
        if window_size > plot_window:
            window_size = plot_window
            print(f"  WARNING: Window size capped at plot_window limit ({plot_window:,} bp)")
        
        # Choose resolution based on whether we had to expand window
        actual_resolution = hic_res_small if use_small_resolution else hic_res
        
        # Calculate region bounds
        half_window = window_size // 2
        region_start = max(0, region_center - half_window)
        region_end = region_center + half_window
        
        # Round to resolution boundaries (using actual resolution for visualization)
        region_start = (region_start // actual_resolution) * actual_resolution
        region_end = ((region_end // actual_resolution) + 1) * actual_resolution
        
        # Create filename
        gene_name = row.get('TargetGene_element1', f'pair_{idx}')
        filename = f"{gene_name}_{chrom}_{region_start}_{region_end}.png"
        save_path = os.path.join(plot_dir, filename)
        
        # Create custom title
        element2_start = row['start_element2']
        title = f"{gene_name}\nE1: {element1_start:,} | E2: {element2_start:,} | TSS: {tss_pos:,}\n{chrom}:{region_start:,}-{region_end:,} (window: {window_size:,}bp, res: {actual_resolution:,}bp)"
        
        print(f"Creating plot {plot_count + 1}/{len(synergy_df)}: {filename}")
        print(f"  Element1: {element1_start:,}, TSS: {tss_pos:,}, Average: {region_center:,}")
        print(f"  E1-TSS distance: {distance_e1_tss:,}, Adaptive: {adaptive_window:,}, Min: {min_window:,}, Final: {window_size:,}")
        print(f"  Using {'small' if use_small_resolution else 'original'} resolution: {actual_resolution:,}")
        
        # Create markers dict for the plot
        markers = {
            'E1': (element1_start, 'blue', 's'),
            'E2': (element2_start, 'green', 's'),
            'TSS': (tss_pos, 'red', 'o')
        }
        
        # Create the Hi-C plot with element markers using chosen resolution
        visualize_hic(hic_file, chrom, region_start, region_end, actual_resolution,
                     save_path=save_path, title=title, markers=markers)
        
        plot_count += 1
        
    except Exception as e:
        print(f"Error creating plot for row {idx}: {e}")
        continue

# Create output file to signal completion
with open(snakemake.output.completed_file, 'w') as f:
    f.write(f"Created {plot_count} Hi-C plots in {plot_dir}\n")

print(f"Successfully created {plot_count} Hi-C visualizations")
