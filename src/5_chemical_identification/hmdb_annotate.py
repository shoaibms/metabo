

######################################################################################
###### Code to Extract All Header (Tag) Names from .zip file and save as .txt ########
######################################################################################

import zipfile
import xml.etree.ElementTree as ET

# Define the path to your .zip file
zip_file = r'C:\Users\ms\Downloads\hmdb_metabolites.zip'

# Create a set to store unique tag names (headers)
unique_headers = set()

# Extract the .xml file from the zip and parse it
with zipfile.ZipFile(zip_file, 'r') as zip_ref:
    xml_file_name = zip_ref.namelist()[0]  # Assuming there's one XML file inside
    print(f"Extracting XML file: {xml_file_name}")

    with zip_ref.open(xml_file_name) as xml_file:
        # Use iterparse to efficiently iterate over the XML structure
        context = ET.iterparse(xml_file, events=("start", "end"))
        
        # Iterate through the XML elements and collect unique tag names
        for event, elem in context:
            if event == "start":
                unique_headers.add(elem.tag)
            elem.clear()  # Clear the element to free memory

# Print all unique headers
print("Unique headers (tags) found in the XML file:")
for header in sorted(unique_headers):
    print(header)

# Optionally, save the headers to a text file
output_file = r'C:\Users\ms\Downloads\hmdb_headers.txt'
with open(output_file, 'w') as f:
    for header in sorted(unique_headers):
        f.write(header + '\n')

print(f"Headers saved to {output_file}")



#################################################################


import zipfile
import xml.etree.ElementTree as ET
import pandas as pd
import os

# File paths
zip_file = r'C:\Users\ms\Downloads\hmdb_metabolites.zip'
output_dir = r'C:/Users/ms/Desktop/data_chem/data/meta/'
output_neg_csv = os.path.join(output_dir, 'negative_matches.csv')
output_pos_csv = os.path.join(output_dir, 'positive_matches.csv')

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Load untargeted data (negative and positive clusters)
neg_cluster_metadata = pd.read_csv(r'C:\Users\ms\Desktop\data_chem\data\meta\neg_cluster_metadata.csv')
pos_cluster_metadata = pd.read_csv(r'C:\Users\ms\Desktop\data_chem\data\meta\pos_cluster_metadata.csv')

# Define a function to safely extract text or return a default value if not present
def safe_extract_text(element, tag, default="Not available"):
    value = element.findtext(tag)
    return value if value is not None and value.strip() != "" else default

# Define the ppm_difference function (this was missing)
def ppm_difference(m1, m2):
    return abs((m1 - m2) / m2) * 1e6

# Extract relevant data from HMDB XML
hmdb_data = []
with zipfile.ZipFile(zip_file, 'r') as zip_ref:
    xml_file_name = zip_ref.namelist()[0]
    with zip_ref.open(xml_file_name) as xml_file:
        context = ET.iterparse(xml_file, events=("start", "end"))
        for event, elem in context:
            if event == "end" and elem.tag == "{http://www.hmdb.ca}metabolite":
                metabolite_data = {
                    'name': safe_extract_text(elem, '{http://www.hmdb.ca}name'),
                    'formula': safe_extract_text(elem, '{http://www.hmdb.ca}chemical_formula'),
                    'exact_mass': safe_extract_text(elem, '{http://www.hmdb.ca}monisotopic_molecular_weight'),
                    'super_class': safe_extract_text(elem, '{http://www.hmdb.ca}super_class'),
                    'sub_class': safe_extract_text(elem, '{http://www.hmdb.ca}sub_class'),
                    'direct_parent': safe_extract_text(elem, '{http://www.hmdb.ca}direct_parent'),
                    'pathways': safe_extract_text(elem, '{http://www.hmdb.ca}pathway'),
                    'biospecimen_locations': safe_extract_text(elem, '{http://www.hmdb.ca}biospecimen_locations'),
                    'tissue_locations': safe_extract_text(elem, '{http://www.hmdb.ca}tissue_locations'),
                    'inchi': safe_extract_text(elem, '{http://www.hmdb.ca}inchi'),
                    'inchikey': safe_extract_text(elem, '{http://www.hmdb.ca}inchikey'),
                    'iupac_name': safe_extract_text(elem, '{http://www.hmdb.ca}iupac_name'),
                    'smiles': safe_extract_text(elem, '{http://www.hmdb.ca}smiles'),
                    'kegg_id': safe_extract_text(elem, '{http://www.hmdb.ca}kegg_id'),
                    'pubchem_compound_id': safe_extract_text(elem, '{http://www.hmdb.ca}pubchem_compound_id'),
                    'kingdom': safe_extract_text(elem, '{http://www.hmdb.ca}kingdom'),
                    'metlin_id': safe_extract_text(elem, '{http://www.hmdb.ca}metlin_id'),
                    'molecular_framework': safe_extract_text(elem, '{http://www.hmdb.ca}molecular_framework'),
                    'synonyms': safe_extract_text(elem, '{http://www.hmdb.ca}synonyms'),
                    'ontology': safe_extract_text(elem, '{http://www.hmdb.ca}ontology'),
                    'parent_id': safe_extract_text(elem, '{http://www.hmdb.ca}parent_id'),
                    'concentration': safe_extract_text(elem, '{http://www.hmdb.ca}concentration'),
                    'predicted_properties': safe_extract_text(elem, '{http://www.hmdb.ca}predicted_properties'),
                    'protein_associations': safe_extract_text(elem, '{http://www.hmdb.ca}protein_associations'),
                    'smpdb_id': safe_extract_text(elem, '{http://www.hmdb.ca}smpdb_id'),
                    'spectrum_id': safe_extract_text(elem, '{http://www.hmdb.ca}spectrum_id'),
                    'substituents': safe_extract_text(elem, '{http://www.hmdb.ca}substituents')
                }

                # Append only if exact_mass is available
                if metabolite_data['exact_mass'] != "Not available":
                    metabolite_data['exact_mass'] = float(metabolite_data['exact_mass'])
                    hmdb_data.append(metabolite_data)
                elem.clear()

# Convert HMDB data to a DataFrame
df_hmdb = pd.DataFrame(hmdb_data)

# Function to calculate m/z and RT ranking
def rank_mz(mz_exp, mz_db):
    return abs(mz_exp - mz_db) / mz_exp

def rank_rt(rt_exp, rt_db, rt_max):
    return abs(rt_exp - rt_db) / rt_max

# Function to match clusters with HMDB data
def match_clusters(cluster_data, mode):
    matches = []
    max_rt = cluster_data['RT'].max()  # Max RT for normalization
    
    for _, row in cluster_data.iterrows():
        mz_value = row['m/z']
        mass_value = row['Mass'] if 'Mass' in row and not pd.isnull(row['Mass']) else None
        rt_value = row['RT']
        
        potential_matches = df_hmdb.copy()

        # If mass is available, prioritize matching based on mass
        if mass_value is not None:
            potential_matches['ppm_diff'] = potential_matches['exact_mass'].apply(lambda x: ppm_difference(mass_value, x))
        else:
            # Fallback to m/z matching if mass is not available
            potential_matches['ppm_diff'] = potential_matches['exact_mass'].apply(lambda x: ppm_difference(mz_value, x))

        # Filter based on 10 ppm tolerance
        valid_matches = potential_matches[potential_matches['ppm_diff'] <= 10]

        if not valid_matches.empty:
            best_match = valid_matches.loc[valid_matches['ppm_diff'].idxmin()]
            
            # Calculate m/z and RT ranks
            mz_rank = rank_mz(mz_value, best_match['exact_mass'])
            rt_rank = rank_rt(rt_value, best_match.get('rt_value', 0), max_rt) if 'rt_value' in best_match else 'Not available'
            
            matches.append({
                f'{mode}_Cluster': row[f'{mode}_Cluster'],
                'Matched Name': best_match['name'],
                'Chemical Formula': best_match['formula'],
                'Chemical Group': best_match['super_class'],
                'Sub Class': best_match['sub_class'],
                'Direct Parent': best_match['direct_parent'],
                'Pathways': best_match['pathways'],
                'Biospecimen Locations': best_match['biospecimen_locations'],
                'Tissue Locations': best_match['tissue_locations'],
                'InChI': best_match['inchi'],
                'InChIKey': best_match['inchikey'],
                'IUPAC Name': best_match['iupac_name'],
                'SMILES': best_match['smiles'],
                'KEGG ID': best_match['kegg_id'],
                'PubChem CID': best_match['pubchem_compound_id'],
                'Kingdom': best_match['kingdom'],
                'Metlin ID': best_match['metlin_id'],
                'Molecular Framework': best_match['molecular_framework'],
                'Synonyms': best_match['synonyms'],
                'Ontology': best_match['ontology'],
                'Parent ID': best_match['parent_id'],
                'Concentration': best_match['concentration'],
                'Predicted Properties': best_match['predicted_properties'],
                'Protein Associations': best_match['protein_associations'],
                'SMPDB ID': best_match['smpdb_id'],
                'Spectrum ID': best_match['spectrum_id'],
                'Substituents': best_match['substituents'],
                'Exact Mass': best_match['exact_mass'],
                'Mass': mass_value if mass_value is not None else 'Not provided',
                'm/z': mz_value,
                'RT': rt_value,
                'ppm_difference': best_match['ppm_diff'],
                'm/z Rank': mz_rank,
                'RT Rank': rt_rank
            })
    
    return pd.DataFrame(matches)

# Perform matching for negative and positive clusters
neg_matches = match_clusters(neg_cluster_metadata, 'N')
pos_matches = match_clusters(pos_cluster_metadata, 'P')

# Save the results to CSV files in the specified directory
neg_matches.to_csv(output_neg_csv, index=False)
pos_matches.to_csv(output_pos_csv, index=False)

print(f"Matching complete. Negative mode matches saved to {output_neg_csv}. Positive mode matches saved to {output_pos_csv}.")