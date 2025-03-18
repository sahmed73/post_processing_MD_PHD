import torch
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdmolops
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader

# Sample antioxidant dataset
data_dict = {
    "SMILES": ["CC1=CC(=O)C=CC1=O",
               "CC(C)CC1=CC(=O)C=CC1=O",
               "Cc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1", # A3
               "CC(C)(C)C1=CC(=CC(=C1O)C(C)(C)C)C=CC(=O)", # A4
               "CC(C)(C)C1=CC(=CC(=C1O)C(C)(C)C)CC2=CC(=C(C(=C2)C(C)(C)C)O)C(C)(C)C" # A5
               ],  
    "Efficiency": [0.85,
                   0.72,
                   0.91, # A3
                   0.87, # A4
                   0.93, # A5
                   ]  # Hypothetical scavenging efficiency
}

df = pd.DataFrame(data_dict)

# Function to convert SMILES to graph
def smiles_to_graph(smiles, target):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # Node features (atomic numbers)
    node_features = torch.tensor([atom.GetAtomicNum() for atom in mol.GetAtoms()], dtype=torch.float)

    # Edge index (connectivity)
    edge_index = []
    for bond in mol.GetBonds():
        edge_index.append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
        edge_index.append([bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()])
    
    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()

    # Convert to PyG Data object
    return Data(x=node_features.view(-1, 1), edge_index=edge_index, y=torch.tensor([target], dtype=torch.float))

# Convert dataset to graphs
graph_data = [smiles_to_graph(smiles, target) for smiles, target in zip(df["SMILES"], df["Efficiency"])]

# DataLoader
loader = DataLoader(graph_data, batch_size=2, shuffle=True)

#%%
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, global_mean_pool

class GCN(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels):
        super().__init__()
        self.conv1 = GCNConv(in_channels, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, hidden_channels)
        self.fc = torch.nn.Linear(hidden_channels, out_channels)  # Final graph-level output

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch  # Batch index needed for pooling
        x = self.conv1(x, edge_index).relu()
        x = self.conv2(x, edge_index).relu()

        # Global pooling: Convert node embeddings to a single graph embedding
        x = global_mean_pool(x, batch)  # (num_graphs, hidden_dim)

        x = self.fc(x)  # Predict efficiency for the whole graph
        return x

#%%# Model setup

model = GCN(in_channels=1, hidden_channels=16, out_channels=1)
optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
loss_fn = torch.nn.MSELoss()  # Regression loss

# Training loop
for epoch in range(100):
    total_loss = 0
    for batch in loader:
        optimizer.zero_grad()
        pred = model(batch).view(-1)
        loss = loss_fn(pred, batch.y)
        loss.backward()
        optimizer.step()
        total_loss += loss.item()

    if epoch % 10 == 0:
        print(f"Epoch {epoch}: Loss = {total_loss:.4f}")

print("Training complete!")

#%%
test_smiles = "C"  # New antioxidant
test_graph = smiles_to_graph(test_smiles, 0.2)  # Placeholder target
test_pred = model(test_graph).item()

print(f"Predicted antioxidant efficiency: {test_pred:.2f}")


