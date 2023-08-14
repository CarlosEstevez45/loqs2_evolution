#!/usr/bin/env python3

#Author: Carlos Estevez-Castro
#Contact: carlosestevez45@ufmg.br / estevezcastro@etu.unistra.fr

#requires pandas v2.0.2 (McKinney et al, https://github.com/pandas-dev/pandas)

import os
import pandas as pd

# Read the TSV table and store normalization factors in a DataFrame
normalization_factors = pd.read_csv('Loqs2_sRNA_larvae_LibsNormValues.tsv', sep='\t', index_col='Lib')

# Define the file patterns to consider
file_patterns = ['miRNAs','siRNAs','piRNAs']

# Process SAM files that match Lib values in the TSV table
for pattern in file_patterns:
    merged_counts = {}

    # Loop through each library
    for lib, norm_factor in normalization_factors.iterrows():
        sam_files = [f for f in os.listdir('.') if f.startswith(f"Loqs2_{lib}_Aae_clusters.FINAL.") and f.endswith(f"{pattern}.sam")]

        # Initialize counts for the current library
        lib_counts = {}

        # Loop through matching files and process them
        for sam_file in sam_files:
            # Read SAM file and extract reference counts
            with open(sam_file, 'r') as f:
                counts = {}
                for line in f:
                    if not line.startswith('@'):
                        ref = line.split('\t')[2]
                        counts[ref] = counts.get(ref, 0) + 1

            # Normalize counts using the corresponding normalization factor
            norm_counts = {ref: count * 1000000 / norm_factor['Norm'] for ref, count in counts.items()}

            # Add normalized counts to the library's count dictionary
            lib_counts.update(norm_counts)

        # Store the library's counts in the merged_counts dictionary
        merged_counts[lib] = lib_counts

    # Convert the merged_counts dictionary to a DataFrame
    merged_counts_df = pd.DataFrame(merged_counts)

    # Add 'Feature' and '1' columns
    merged_counts_df.insert(0, 'Feature', merged_counts_df.index)

    # Replace NaN values with 0
    merged_counts_df.fillna(0, inplace=True)

    # Save merged counts to a TSV file for the current file pattern
    output_filename = f'merged_NormalizedCounts_{pattern}.tsv'
    merged_counts_df.to_csv(output_filename, sep='\t', index=False)
    print(f"Merged normalized counts for {pattern} saved to {output_filename}")

    import subprocess
    # specify the path to your R script
    r_script_path = "./HMofNormCounts.R" # replace with actual path
    # define your filtering values
    filtering_values = [5,50,100, 500]

    # run the R script with each filtering value
    for filtering in filtering_values:
        # call the R script with the current filtering value
        retcode = subprocess.call(["Rscript", r_script_path, output_filename, str(filtering)])
        if retcode != 0:
            print(f"R script ended with error code {retcode}")
