from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import pandas as pd

# Sample data: SMILES strings
data = {
    'smiles': ['CCO', 'CC(=O)O', 'CCOCC', 'CCN', 'CCOCCN', 'CCCC', 'CCOC'],
    'name': ['ethanol', 'acetic acid', 'ethyl ether', 'ethylamine', 'diethylamine', 'butane', 'methoxyethane']
}
df = pd.DataFrame(data)

# Function to generate Morgan fingerprint
def generate_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

# Generate fingerprints for all molecules in the dataset
df['fingerprint'] = df['smiles'].apply(generate_fingerprint)

# Function to calculate similarity
def calculate_similarity(fp1, fp2):
    return DataStructs.FingerprintSimilarity(fp1, fp2)

# Example: Find the most similar compound to 'ethyl ether'
target_smiles = 'CCOCC'
target_fp = generate_fingerprint(target_smiles)

# Calculate similarity to all other compounds
df['similarity'] = df['fingerprint'].apply(lambda x: calculate_similarity(target_fp, x))

# Sort by similarity
df_sorted = df.sort_values(by='similarity', ascending=False)

# Display results
print(df_sorted[['name', 'smiles', 'similarity']])
