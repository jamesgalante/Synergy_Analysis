## Get Hi-C interaction frequencies for all E-G pairs on the same chromosome
import hicstraw
import pandas as pd

# Define functions ---------------------------------------------------------------------------------

# look up contact frequency for one pair using full bin ranges
def get_int_freq_pair(bin1_start, bin1_end, bin2_start, bin2_end, matrix_chr):
    int_freq = matrix_chr.getRecords(int(bin1_start), int(bin1_end), int(bin2_start), int(bin2_end))
    if len(int_freq) > 0:
        int_freq = int_freq[0].counts
    else:
        int_freq = float("nan")
    return(int_freq)

# look up all contact frequencies for E-G pairs on one chromosome
def get_eg_int_freqs_chr(etp, hic, chrom, hic_res):
    
    if len(etp) == 0:
        return etp
    
    # Initialize all result columns with NaN
    contact_columns = ['e1_tss_obs_scale', 'e2_tss_obs_scale', 'e1_e2_obs_scale']
    
    for col in contact_columns:
        etp[col] = float('nan')
    
    # Process one matrix type at a time to save memory
    print("    Processing observed + SCALE")
    matrix = hic.getMatrixZoomData(chrom, chrom, "observed", "SCALE", "BP", hic_res)
    
    for idx, row in etp.iterrows():
        etp.at[idx, 'e1_tss_obs_scale'] = get_int_freq_pair(
            row['enh_hic_start_element1'], row['enh_hic_end_element1'],  
            row['tss_hic_start'], row['tss_hic_end'], matrix)
        etp.at[idx, 'e2_tss_obs_scale'] = get_int_freq_pair(
            row['enh_hic_start_element2'], row['enh_hic_end_element2'], 
            row['tss_hic_start'], row['tss_hic_end'], matrix)
        etp.at[idx, 'e1_e2_obs_scale'] = get_int_freq_pair(
            row['enh_hic_start_element1'], row['enh_hic_end_element1'], 
            row['enh_hic_start_element2'], row['enh_hic_end_element2'], matrix)
    del matrix
    
    return etp

# get contact frequencies for all E-G pairs in a ETP table (single chromosome)
def get_eg_int_freqs(etp, hic, hic_res):
    
    # Element 1: calculate hic_bin for each enhancer based on enhancer center
    etp['enh_center_element1'] = (etp['start_element1'] + etp['end_element1']) / 2
    etp['enh_hic_start_element1'] = etp['enh_center_element1'] // hic_res * hic_res
    etp['enh_hic_end_element1'] = etp['enh_hic_start_element1'] + hic_res
    
    # Element 2: calculate hic_bin for each enhancer based on enhancer center
    etp['enh_center_element2'] = (etp['start_element2'] + etp['end_element2']) / 2
    etp['enh_hic_start_element2'] = etp['enh_center_element2'] // hic_res * hic_res
    etp['enh_hic_end_element2'] = etp['enh_hic_start_element2'] + hic_res
    
    # calculate hic_bin for each gene TSS (one common TSS for each element - just using element1)
    etp['tss_hic_start'] = etp['TargetGeneTSS_element1'] // hic_res * hic_res
    etp['tss_hic_end'] = etp['tss_hic_start'] + hic_res
    
    # Process this single chromosome
    chrom = etp['#chr_element1'].iloc[0]
    print(f"Getting interaction frequencies for chromosome {chrom}")
    etp_hic = get_eg_int_freqs_chr(etp, hic, chrom, hic_res)
    
    return(etp_hic)

# Get contact frequencies for E-G pairs ------------------------------------------------------------

# load ETP table
etp = pd.read_csv(snakemake.input.paired_predictions, low_memory = False, sep = "\t")

# Filter for current chromosome
etp = etp[etp['#chr_element1'] == snakemake.wildcards.chrom].copy()

# open connection to .hic file
hic = hicstraw.HiCFile(snakemake.input.hic)

# Print available resolutions for debugging
available_resolutions = hic.getResolutions()
requested_resolution = snakemake.params.hic_res

print(f"Available resolutions in Hi-C file: {available_resolutions}")
print(f"Requested resolution: {requested_resolution}")

if requested_resolution not in available_resolutions:
    print(f"ERROR: Requested resolution {requested_resolution} not available!")
    print(f"Available options: {available_resolutions}")
    raise ValueError(f"Resolution {requested_resolution} not found in Hi-C file")

print(f"Starting analysis with {len(etp)} enhancer pairs...")

# get interaction frequencies for all E-G pairs in ETP table and write output to file
etp_hic = get_eg_int_freqs(etp, hic, hic_res = requested_resolution)

# Save output - drop intermediate columns
columns_to_drop = [
    'enh_center_element1', 'enh_hic_start_element1', 'enh_hic_end_element1',
    'enh_center_element2', 'enh_hic_start_element2', 'enh_hic_end_element2', 
    'tss_hic_start', 'tss_hic_end'
]

# Only drop columns that actually exist
columns_to_drop = [col for col in columns_to_drop if col in etp_hic.columns]
if columns_to_drop:
    etp_hic = etp_hic.drop(columns=columns_to_drop)

print(f"Saving results with {len(etp_hic)} rows and contact columns: {[col for col in etp_hic.columns if any(x in col for x in ['e1_', 'e2_'])]}")
etp_hic.to_csv(snakemake.output.paired_predictions_with_hic_per_chrom, index = False, sep = "\t")
