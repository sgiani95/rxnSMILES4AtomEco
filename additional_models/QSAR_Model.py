from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.ensemble import RandomForestRegressor
import pandas as pd

# Sample data: SMILES strings and their corresponding activity/property values
data = {
    'smiles': ['CCO', 'CC(=O)O', 'CCOCC', 'CCN'],
    'activity': [1.2, 2.3, 0.8, 3.1]
}
df = pd.DataFrame(data)

# Calculate molecular descriptors
def calculate_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return [Descriptors.MolWt(mol), Descriptors.TPSA(mol), Descriptors.MolLogP(mol)]

df['descriptors'] = df['smiles'].apply(calculate_descriptors)
df = df.dropna()  # Drop rows where descriptor calculation failed

# Prepare data for model training
X = pd.DataFrame(df['descriptors'].tolist())
y = df['activity']

# Train a RandomForest model
model = RandomForestRegressor()
model.fit(X, y)

# Predict activity for a new compound
new_smiles = 'CCOCCN'
new_descriptors = calculate_descriptors(new_smiles)
predicted_activity = model.predict([new_descriptors])
print(f'Predicted activity for {new_smiles}: {predicted_activity[0]}')
