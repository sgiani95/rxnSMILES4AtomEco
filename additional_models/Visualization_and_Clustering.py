import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, DataStructs
from sklearn.cluster import KMeans
import pandas as pd
import matplotlib.pyplot as plt

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

# Convert RDKit fingerprints to numpy arrays for clustering
def fingerprint_to_numpy(fp):
    arr = np.zeros((1,), dtype=np.int32)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

fingerprint_matrix = np.array([fingerprint_to_numpy(fp) for fp in df['fingerprint']])

# Perform KMeans clustering
n_clusters = 3
kmeans = KMeans(n_clusters=n_clusters, random_state=42).fit(fingerprint_matrix)
df['cluster'] = kmeans.labels_

# Visualize clusters
for cluster in range(n_clusters):
    cluster_data = df[df['cluster'] == cluster]
    mols = [Chem.MolFromSmiles(smiles) for smiles in cluster_data['smiles']]
    img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(200, 200))
    plt.figure(figsize=(10, 10))
    plt.title(f'Cluster {cluster + 1}')
    plt.imshow(img)
    plt.axis('off')
    plt.show()
