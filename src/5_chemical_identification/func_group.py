import os
import pandas as pd
from rdkit import Chem
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

# Comprehensive SMARTS patterns for functional groups in plant metabolites
functional_groups = {
    'Hydroxyl': '[OX2H]',                           # Alcohols
    'Phenol': 'c1ccc(O)cc1',                        # Phenols
    'Carboxyl': 'C(=O)[OH]',                        # Carboxylic acids
    'Aldehyde': '[CX3H1](=O)',                      # Aldehydes
    'Ketone': 'C(=O)[C]',                           # Ketones
    'Ester': 'C(=O)O[C]',                           # Esters
    'Ether': '[OD2]([#6])[#6]',                     # Ethers
    'Amine': '[NX3;H2,H1;!$(NC=O)]',                # Amines
    'Amide': 'C(=O)N',                              # Amides
    'Alkene': 'C=C',                                # Alkenes
    'Alkyne': 'C#C',                                # Alkynes
    'Aromatic': 'c1ccccc1',                         # Aromatic rings (Benzene)
    'Terpene': '[CX2H]C=C',                         # Terpenes (isoprenoids)
    'Alkaloid': '[nH]1cccc1',                       # Alkaloids (heterocyclic amines)
    'Flavonoid': 'c1cc(O)c2c(c1)cc(O)cc2',          # Flavonoids (polyphenols)
    'Halogen': '[F,Cl,Br,I]',                       # Halogenated compounds
    'Nitro': '[$([NX3](=O)=O)]',                    # Nitro compounds
    'Sulfonyl': 'S(=O)(=O)[O]',                     # Sulfonyl group
    'Thiol': '[SX2H]',                              # Thiols
    'Phosphate': 'P(=O)(O)(O)',                     # Phosphate group
    'Glucoside': '[C,CX4]O[C,CX4]',                 # Glycosides (Glucosides)
    'Isocyanide': '[NX1]#[CX1]',                    # Isocyanides
    'Epoxide': 'C1OC1',                             # Epoxides (cyclic ethers)
    'Quinone': 'O=C1C=CC(=O)C=C1',                  # Quinones
    'Isothiocyanate': 'N=C=S',                      # Isothiocyanates
    'Coumarin': 'O=C1OCc2ccccc12',                  # Coumarins
    'Pyrrole': 'C1=CC=CN1',                         # Pyrroles
    'Furan': 'C1=COC=C1',                           # Furans
    'Pyrazole': 'C1=CC=NN1',                        # Pyrazoles
    'Pyridine': 'C1=CC=NC=C1',                      # Pyridines
}

# Function to extract functional groups from a SMILES string
def get_functional_groups(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "Invalid SMILES"
        
        detected_groups = []
        for group_name, smarts in functional_groups.items():
            if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
                detected_groups.append(group_name)
        
        return ", ".join(detected_groups) if detected_groups else "No functional groups"
    except Exception as e:
        return f"Error: {str(e)}"

# Function to process individual files
def process_file(input_file, output_file):
    df = pd.read_csv(input_file)
    with ProcessPoolExecutor() as executor:
        df['Functional_Groups'] = list(tqdm(
            executor.map(get_functional_groups, df['SMILES']),
            total=len(df),
            desc=f"Processing {os.path.basename(input_file)}"
        ))
    df.to_csv(output_file, index=False)
    print(f"File saved: {output_file}")

# Main function to orchestrate the processing of input files
def main():
    base_dir = r"C:\Users\USER\Desktop\data_chem_3_10\output\results"
    input_files = {
        "negative": os.path.join(base_dir, "rdkit", "negative_matches_with_fingerprints_with_scaled_confidence.csv"),
        "positive": os.path.join(base_dir, "rdkit", "positive_matches_with_fingerprints_with_scaled_confidence.csv")
    }
    output_dir = os.path.join(base_dir, "group", "functional_group")
    os.makedirs(output_dir, exist_ok=True)

    for file_type, input_file in input_files.items():
        output_file = os.path.join(output_dir, f"{file_type}_matches_with_functional_groups.csv")
        process_file(input_file, output_file)

# Entry point for the script
if __name__ == "__main__":
    main()