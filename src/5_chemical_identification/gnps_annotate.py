
#############################################
########### GNPS Search #####################
#############################################

import pandas as pd
import requests
import os

# File paths
neg_cluster_metadata = r'C:\Users\ms\Desktop\data_chem\data\meta\neg_cluster metadata.csv'
pos_cluster_metadata = r'C:\Users\ms\Desktop\data_chem\data\meta\pos_cluster metadata.csv'
negative_matches = r'C:\Users\ms\Desktop\data_chem\data\meta\negative_matches.csv'
positive_matches = r'C:\Users\ms\Desktop\data_chem\data\meta\positive_matches.csv'

output_neg_csv = r'C:\Users\ms\Desktop\data_chem\data\meta\negative_matches_enriched_gnps.csv'
output_pos_csv = r'C:\Users\ms\Desktop\data_chem\data\meta\positive_matches_enriched_gnps.csv'

# Your GNPS API key
gnps_api_key = 'f0d5aa78-5771-40b6-932c-5b6a25f2a3fc'

# Function to query GNPS using InChIKey, SMILES, or PubChem CID
def query_gnps(identifier, identifier_type):
    url = f"https://gnps.ucsd.edu/gnpsapi/metabolites/{identifier_type}/{identifier}"
    headers = {"Authorization": f"Bearer {gnps_api_key}"}
    response = requests.get(url, headers=headers)
    
    if response.status_code == 200:
        try:
            data = response.json()
            chemical_group = data.get('chemical_classification', {}).get('Class', 'Not available')
            pathways = data.get('pathways', 'Not available')
            biospecimen = data.get('biospecimen_locations', 'Not available')
            synonyms = data.get('synonyms', ['Not available'])[0]
            
            return {
                'Chemical Group': chemical_group,
                'Pathways': pathways,
                'Biospecimen Locations': biospecimen,
                'Synonyms': synonyms
            }
        except (IndexError, KeyError):
            return {
                'Chemical Group': 'Not available',
                'Pathways': 'Not available',
                'Biospecimen Locations': 'Not available',
                'Synonyms': 'Not available'
            }
    else:
        return {
            'Chemical Group': 'Not available',
            'Pathways': 'Not available',
            'Biospecimen Locations': 'Not available',
            'Synonyms': 'Not available'
        }

# Function to enrich the data with GNPS API using InChIKey, SMILES, or PubChem CID
def enrich_with_gnps(df):
    enriched_data = []
    for _, row in df.iterrows():
        if pd.notna(row['InChIKey']) and row['InChIKey'] != 'Not available':
            gnps_info = query_gnps(row['InChIKey'], 'inchikey')
        elif pd.notna(row['SMILES']) and row['SMILES'] != 'Not available':
            gnps_info = query_gnps(row['SMILES'], 'smiles')
        elif pd.notna(row['PubChem CID']) and row['PubChem CID'] != 'Not available':
            gnps_info = query_gnps(row['PubChem CID'], 'cid')
        else:
            gnps_info = {
                'Chemical Group': 'Not available',
                'Pathways': 'Not available',
                'Biospecimen Locations': 'Not available',
                'Synonyms': 'Not available'
            }
        
        enriched_data.append(gnps_info)
    
    # Convert the enriched data to DataFrame
    enriched_df = pd.DataFrame(enriched_data)
    enriched_full_df = pd.concat([df, enriched_df], axis=1)
    return enriched_full_df

# Load the negative and positive match files
neg_matches_df = pd.read_csv(negative_matches)
pos_matches_df = pd.read_csv(positive_matches)

# Enrich the data with GNPS API
enriched_neg_matches = enrich_with_gnps(neg_matches_df)
enriched_pos_matches = enrich_with_gnps(pos_matches_df)

# Save the enriched data as CSV
enriched_neg_matches.to_csv(output_neg_csv, index=False)
enriched_pos_matches.to_csv(output_pos_csv, index=False)

print(f"Enriched negative matches saved to: {output_neg_csv}")
print(f"Enriched positive matches saved to: {output_pos_csv}")



#########################################################################
##### incorporates tissue type information in the above file ############
#########################################################################



import pandas as pd

# Define file paths
negative_file_path = r'C:\\Users\\ms\\Desktop\\data_chem_3_10\\data\\meta\\negative_matches_enriched_gnps.csv'
positive_file_path = r'C:\\Users\\ms\\Desktop\\data_chem_3_10\\data\\meta\\positive_matches_enriched_gnps.csv'
cluster_meta_file_path = r'C:\\Users\\ms\\Desktop\\data_chem_3_10\\data\\old_1\\Cluster_meta_tissue_allocation.csv'
negative_output_file_path = r'C:\\Users\\ms\\Desktop\\data_chem_3_10\\data\\meta\\negative_matches_enriched_gnps_tissue_type.csv'
positive_output_file_path = r'C:\\Users\\ms\\Desktop\\data_chem_3_10\\data\\meta\\positive_matches_enriched_gnps_tissue_type.csv'

# Load the CSV and metadata files
negative_matches = pd.read_csv(negative_file_path)
positive_matches = pd.read_csv(positive_file_path)
cluster_meta = pd.read_csv(cluster_meta_file_path)

# Extract clusters related to root and leaf for both negative and positive modes
n_clusters_root = cluster_meta['N_Cluster_root'].dropna().unique()
n_clusters_leaf = cluster_meta['N_Cluster_leaf'].dropna().unique()
p_clusters_root = cluster_meta['P_Cluster_root'].dropna().unique()
p_clusters_leaf = cluster_meta['P_Cluster_leaf'].dropna().unique()

# Add root and leaf columns to negative matches
negative_matches['root'] = negative_matches['N_Cluster'].apply(lambda x: 'yes' if x in n_clusters_root else 'no')
negative_matches['leaf'] = negative_matches['N_Cluster'].apply(lambda x: 'yes' if x in n_clusters_leaf else 'no')

# Add root and leaf columns to positive matches
positive_matches['root'] = positive_matches['P_Cluster'].apply(lambda x: 'yes' if x in p_clusters_root else 'no')
positive_matches['leaf'] = positive_matches['P_Cluster'].apply(lambda x: 'yes' if x in p_clusters_leaf else 'no')

# Save the updated files
negative_matches.to_csv(negative_output_file_path, index=False)
positive_matches.to_csv(positive_output_file_path, index=False)

print("Files updated and saved successfully.")






#####################################################################################
##### incorporates VIP_mann_whitney_bonferroni_fdr with the above files #############
############ updated code with corss validation #####################################
#####################################################################################

import pandas as pd

# Define file paths
negative_output_file_path = r'C:\\Users\\ms\\Desktop\\data_chem_3_10\\data\\meta\\negative_matches_enriched_gnps_tissue_type.csv'
positive_output_file_path = r'C:\\Users\\ms\\Desktop\\data_chem_3_10\\data\\meta\\positive_matches_enriched_gnps_tissue_type.csv'
vip_mann_whitney_file_path = r'C:\Users\ms\Desktop\data_chem_3_10\output\results\vip_bonferroni\VIP_mann_whitney_bonferroni_fdr_combine.csv'

# Load the CSV files and the VIP results file
negative_matches = pd.read_csv(negative_output_file_path)
positive_matches = pd.read_csv(positive_output_file_path)
vip_mann_whitney_data = pd.read_csv(vip_mann_whitney_file_path)

# Merge the VIP data with negative_matches and positive_matches based on N_Cluster and P_Cluster
# Negative mode
negative_matches = pd.merge(negative_matches, vip_mann_whitney_data, 
                            left_on='N_Cluster', right_on='Metabolite', how='left')

# Positive mode
positive_matches = pd.merge(positive_matches, vip_mann_whitney_data, 
                            left_on='P_Cluster', right_on='Metabolite', how='left')

# Save the updated files
negative_output_file_path_updated = r'C:\\Users\\ms\\Desktop\\data_chem_3_10\\data\\meta\\negative_matches_enriched_gnps_final.csv'
positive_output_file_path_updated = r'C:\\Users\\ms\\Desktop\\data_chem_3_10\\data\\meta\\positive_matches_enriched_gnps_final.csv'

negative_matches.to_csv(negative_output_file_path_updated, index=False)
positive_matches.to_csv(positive_output_file_path_updated, index=False)

print("Files updated and saved successfully.")