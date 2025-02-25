

from scipy.io import loadmat
from scipy.io import savemat
import numpy as np


data = loadmat("../data/Reconstructions_data/Roseburia_intestinalis_L1_82.mat")


reactions = data['model'][0][0]['rxns']
reaction_ids = [str(r[0]) for r in reactions]
duplicates = [r for r in set(reaction_ids) if reaction_ids.count(r) > 1]
print("Duplicate reactions in model:" , duplicates)


metabolites = data['model'][0][0]['mets']
metabolites_ids = [str(r[0]) for r in metabolites]
duplicates = [r for r in set(metabolites_ids) if metabolites_ids.count(r) > 1]
print("Duplicate metabolites in model:" , duplicates)


genes = data['model'][0][0]['genes']
genes_ids = [r[0].strip("[]'") if isinstance(r[0], str) else str(r[0]).strip("[]'") for r in genes]
duplicates = [r for r in set(genes_ids) if genes_ids.count(r) > 1]
print("Duplicate genes in model:" , duplicates)


fixed_genes_ids = []
seen = {}

for r in genes_ids:
    if r in seen:
        # add a suffix to make the ID unique
        unique_id = f"{r}_{seen[r]}"
        seen[r] += 1
        fixed_genes_ids.append(unique_id)
    else:
        seen[r] = 1
        fixed_genes_ids.append(r)


# verify the fixed IDs
#print(fixed_genes_ids)


# modify the names in the model
data['model'][0][0]['genes'] = np.array([[r] for r in fixed_genes_ids], dtype=object)


# save the model
#savemat("../data/Reconstructions_data/Roseburia_intestinalis_L1_82_fixed.mat", data)
