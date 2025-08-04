"""
Chemical Structure Classification using SMARTS Patterns.

This script classifies chemical compounds from a CSV file based on their SMILES strings.
It uses SMARTS patterns to identify the chemical class of each compound, calculates
a confidence score for the classification, and assigns the most likely class based
on a predefined priority order.

The script requires RDKit and pandas to be installed.

Note to users: Please update the file paths in the CONFIGURATION section to match
your local environment before running the script.
"""

import os
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd
from rdkit import Chem

# --- CONFIGURATION ---

# NOTE: The user requested to keep these hardcoded paths.
# You will need to update these paths to match your local system setup.
NEGATIVE_FILE_PATH = r'C:\Users\ms\Desktop\data_chem_3_10\output\results\rdkit\positive_matches_with_fingerprints.csv'
OUTPUT_DIR = r'C:\Users\ms\Desktop\data_chem_3_10\output\results\rdkit'

# Priority order for classification, based on abundance and detectability in LCMS.
# This list is used to resolve ties when a compound matches multiple classes
# with the same score.
PRIORITY_ORDER = [
    "Lipids", "Carbohydrates", "Phenolics", "Amino Acids", "Glycosides",
    "Terpenoids", "Proteins", "Alkaloids", "Oxylipins", "Phytochelatins",
    "Volatile Organic Compounds (VOCs)"
]

# Manual weights for specific biomolecule groups. Higher weights give more priority
# during the scoring process.
GROUP_WEIGHTS = {
    "Lipids": 1.2,
    "Carbohydrates": 1.1,
}

# SMARTS patterns for general biomolecule classes.
BIOMOLECULE_SMARTS = {
    "Carbohydrates": [
        "[C@H]([OH])[C@H](O)[C@H](O)[C@H](O)[C@H](O)[CH2OH]",
        "[OX2H][CX4H]1[OX2][CX4H][CX4H][CX4H][CX4H]1[OX2H]",
        "[C@H]([OH])[C@H](O)[C@H](O)[C@H](O)[C@H](O)[OX2H]"
    ],
    "Lipids": [
        "[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2]C(=O)O",
        "C(C(CO)O)O",
        "[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][C](=O)[O]",
        "[C](=O)[O]"
    ],
    "Proteins": [
        "[NX3][CX4H]([*])[CX3](=[OX1])[NX3H]",
        "[NX3][CX4H]([*])[CX3](=[OX1])[OX2H]",
        "[#6](-[#7])-[#6](=[#8])-[#7]",
        "[N][C][C](=O)[N]"
    ],
    "Amino Acids": [
        "[NX3H2][CX4H]([*])[CX3](=[OX1])[OX2H]",
        "[#7][#6](-[*])[#6](=[#8])[#8]",
        "[NH2][CH](*)[C](=O)[OH]",
        "N[C@@H](C(=O)O)"
    ],
    "Phenolics": [
        "[OX2H][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1",
        "[OX2H][c]1[c][c][c][c][c]1",
        "[OX2H][c]1[c][c]([OX2H])[c][c][c]1"
    ],
    "Terpenoids": [
        "[CH3][C](=[CH2])[CH2][CH2][C](=[CH2])[CH3]",
        "[C](=[CH2])[C](=[CH2])[C](=[CH2])[C](=[CH2])",
        "[C@H]1[C@H]([C@H]([C@H](C1)C)[CH2])[CH2]"
    ],
    "Alkaloids": [
        "[$([NX3H2]),$([NX3H](C)(C)),$([NX4H2+](C)(C)(C))]",
        "[#7]",
        "[#6]~[#7]~[#6]~[#7]",
        "C1CC[NH+](C1)C2=CC=CC=C2"
    ],
    "Glycosides": [
        "[C,O][C,O][C,O][C,O][C,O][OX2R0,OX2R1]",
        "[#6]1[#8][#6]([#8])[#6]([#8])[#6]([#8])[#6]1[#8]",
        "[C@H]1([C@H]([C@H](O[C@H]1O)[OH])[OH])[OH]"
    ],
    "Volatile Organic Compounds (VOCs)": [
        "[#6][#6][#6]",
        "[C,H][C,H][C,H]",
        "[#6]~[#6]~[#8,#7]",
        "[CX3H1](=O)[CX3H2]"
    ],
    "Phytochelatins": [
        "[NX3][CX4H]([CH2]S)[CX3](=[OX1])[NX3H]",
        "[NX3][CX4H]([*])[CX3](=[OX1])[NX3H]",
        "[#7][#6](-[#6]-[#16])[#6](=[#8])[#7]"
    ],
    "Oxylipins": [
        "[CH2][CH2][CH2]C=O",
        "[CX4H2][CX3](=[OX1])[CX4H2][CX3H]=[CX3H]",
        "[#6]~[#6]~[#6]~[#6]~[#6](=[#8])[#8]"
    ]
}

# SMARTS patterns for specific, reviewed compounds.
REVIEWED_SMARTS = {
    "Abscisic acid": "CC(=C(C(=O)O)C1=CC(=O)C(=CC1=O)C(C)C)C",
    "Brassinolids": "[#6H3]-[#6H](-[#6H3])-[#6@H](-[#6H3])-[#6@H](-[#6@@H](-[#6@@H](-[#6H3])-[#6@H]1-[#6H2]-[#6H2]-[#6@H]2-[#6@@H]3-[#6H2]-[#8]-[#6](-[#6@H]4-[#6H2]-[#6@@H](-[#6@@H](-[#6H2]-[#6@]-4(-[#6H3])-[#6@H]-3-[#6H2]-[#6H2]-[#6@]-1-2-[#6H3])-[#8H])-[#8H])=[#8])-[#8H])-[#8H]",
    "Ethylene": "C=C",
    "Salicylic Acid": "OC(=O)c1ccccc1O",
    "Jasmonates": "CC(=O)CCC=CC=CC(C)C(=O)O",
    "Auxins": "C1=CC=C2C(=C1)C(=O)NC2=O",
    "Gibberellins": "C1=CC2=C(C=C1)C3C(C2)C4C(C3)C5C(C4)C(C5)C(=O)O",
    "Cytokinins": "C1=NC2=C(N1)C(=NC=N2)N",
    "Strigolactones": "CC1=CC2=C(C=C1)C(=O)OC2=O",
    "Glycerolipids": "C(C(CO)O)O",
    "Phospholipids": "C(C(COP(=O)(O)O)O)O",
    "Galactolipids": "C1C(C(C(C(C(O1)O)O)O)O)CO",
    "Sphingolipids": "NC(CO)C(O)CCCCCCCCCCCCC=C",
    "Superoxide anion": "[O-][O]",
    "Hydrogen peroxide": "OO",
    "Hydroxyl radical": "[OH]",
    "Singlet oxygen": "O=[O]",
    "Fructose": "C(C(C(C(C=O)O)O)O)O",
    "Glucose": "C(C1C(C(C(C(O1)O)O)O)O)O",
    "Proline": "C1CC(NC1)C(=O)O",
    "Lignin": [
        "[C]1=C[C]=C[C]=C1",
        "[C]1=C([O])C=C([C])C=C1",
        "[C]1=C([O])C([O])=C([C])C([O])=C1",
        "[C]1=C([O])C=C([C])C=C1"
    ],
    "Malic acid": "OC(CC(O)=O)C(O)=O"
}

# --- Functions ---

def load_data(file_path: str) -> pd.DataFrame:
    """
    Loads data from a CSV file into a pandas DataFrame.

    Args:
        file_path: The path to the CSV file.

    Returns:
        A pandas DataFrame with the loaded data.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Input file not found at: {file_path}")
    return pd.read_csv(file_path)


def get_smarts_score(smiles: str, smarts_pattern: Union[str, List[str]]) -> Tuple[int, int, int]:
    """
    Calculates a score for a SMILES string based on SMARTS patterns.

    Args:
        smiles: The SMILES string of the molecule.
        smarts_pattern: A single SMARTS string or a list of SMARTS strings.

    Returns:
        A tuple containing the best complexity score, number of atoms matched,
        and the confidence score from the best matching SMARTS pattern.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return 0, 0, 0

    if isinstance(smarts_pattern, str):
        smarts_pattern = [smarts_pattern]

    best_complexity = 0
    best_atoms_matched = 0
    best_confidence = 0

    for smarts in smarts_pattern:
        substruct = Chem.MolFromSmarts(smarts)
        if substruct and mol.HasSubstructMatch(substruct):
            match = mol.GetSubstructMatch(substruct)
            num_atoms_matched = len(match)
            complexity_score = len(smarts)
            confidence = complexity_score * num_atoms_matched

            if confidence > best_confidence:
                best_complexity = complexity_score
                best_atoms_matched = num_atoms_matched
                best_confidence = confidence

    return best_complexity, best_atoms_matched, best_confidence


def best_match_smarts(
    smiles: str,
    all_smarts: Dict[str, Union[str, List[str]]],
    priority_order: List[str],
    group_weights: Optional[Dict[str, float]] = None
) -> Tuple[str, float]:
    """
    Determines the best chemical classification for a SMILES string.

    Args:
        smiles: The SMILES string of the molecule.
        all_smarts: A dictionary mapping group names to SMARTS patterns.
        priority_order: A list of group names for tie-breaking.
        group_weights: A dictionary of weights for specific groups.

    Returns:
        A tuple containing the best classification group and its scaled confidence score.
    """
    if not smiles or pd.isna(smiles):
        return "Unclassified", 0.0

    group_weights = group_weights or {}
    all_matches = []
    all_confidences = []

    for group, smarts_pattern in all_smarts.items():
        complexity, atoms_matched, confidence = get_smarts_score(smiles, smarts_pattern)
        all_confidences.append(confidence)

        if atoms_matched > 0:
            weight = group_weights.get(group, 1.0)
            score = (complexity * weight) + atoms_matched
            all_matches.append({"group": group, "score": score, "confidence": confidence})

    if not all_matches:
        return "Unclassified", 0.0

    best_score = max(match['score'] for match in all_matches)
    top_matches = [match for match in all_matches if match['score'] == best_score]

    if len(top_matches) > 1:
        top_matches.sort(key=lambda m: priority_order.index(m['group']) if m['group'] in priority_order else float('inf'))

    best_match = top_matches[0]
    best_group = best_match['group']
    best_confidence = best_match['confidence']

    min_conf = min(all_confidences)
    max_conf = max(all_confidences)

    scaled_confidence = 0.0
    if max_conf > min_conf:
        scaled_confidence = (best_confidence - min_conf) / (max_conf - min_conf)
    elif max_conf > 0:
        scaled_confidence = 1.0

    return best_group, scaled_confidence


def classify_molecules(df: pd.DataFrame, all_smarts: Dict, priority: List, weights: Dict) -> pd.DataFrame:
    """
    Applies the classification to a DataFrame of molecules.

    Args:
        df: DataFrame containing a 'SMILES' column.
        all_smarts: Combined dictionary of SMARTS patterns.
        priority: The priority order for tie-breaking.
        weights: Group-specific weights.

    Returns:
        The DataFrame with added 'Best_Classification' and 'Scaled_Confidence_Score' columns.
    """
    results = df['SMILES'].apply(lambda x: best_match_smarts(x, all_smarts, priority, weights))
    df[['Best_Classification', 'Scaled_Confidence_Score']] = pd.DataFrame(results.tolist(), index=df.index)
    return df


def save_results(df: pd.DataFrame, output_dir: str, filename: str) -> str:
    """
    Saves the DataFrame to a CSV file.

    Args:
        df: The DataFrame to save.
        output_dir: The directory to save the file in.
        filename: The name of the output file.

    Returns:
        The full path to the saved file.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, filename)
    df.to_csv(output_path, index=False)
    return output_path


def main():
    """
    Main function to run the chemical structure classification workflow.
    """
    print("Starting chemical classification...")
    try:
        df = load_data(NEGATIVE_FILE_PATH)
        
        all_smarts_patterns = {**BIOMOLECULE_SMARTS, **REVIEWED_SMARTS}
        
        classified_df = classify_molecules(df, all_smarts_patterns, PRIORITY_ORDER, GROUP_WEIGHTS)
        
        output_filename = 'positive_matches_with_fingerprints_with_scaled_confidence.csv'
        saved_path = save_results(classified_df, OUTPUT_DIR, output_filename)
        
        print(f"Classification complete. File saved at: {saved_path}")

    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Please ensure the file paths in the script are correct.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    main()
