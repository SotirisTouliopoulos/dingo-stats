

from statsmodels.stats.multitest import multipletests
from scipy import stats
import numpy as np
        

def significantly_altered_reactions(conditions=[], selected_comparisons = [(0, 1)], cobra_model=None, p_value_cutoff=0.05):
    cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]
        
    # if provided sampling dataset has reactions as cols ==> transpose
    for condition in conditions:
        if condition.shape[1] == len(cobra_reactions_str):
            condition = condition.T
    
    # Store p-values for each row for KS test
    p_values_ks = {row: [] for row in range(len(cobra_reactions_str)) }
    # Store ks-values for each row for KS test
    values_ks = {row: [] for row in range(len(cobra_reactions_str)) }

    for row in range(len(cobra_reactions_str)):
        for i, j in selected_comparisons:
            
            ks, p_ks = stats.ks_2samp(conditions[i][row], conditions[j][row], alternative='two-sided')

            p_values_ks[row].append(p_ks)
            values_ks[row].append(ks)
    
    p_values_ks_copy = list(p_values_ks.values())
    flat_p_values_ks = np.array(p_values_ks_copy).flatten()

    # Apply FDR correction
    _, corrected_p_values_ks, _, _ = multipletests(flat_p_values_ks, method='fdr_bh')

    # Reshape corrected p-values back to the original matrix shape
    corrected_p_values_ks = corrected_p_values_ks.reshape(np.array(p_values_ks_copy).shape)
    
    significant_diff_indices = np.where(corrected_p_values_ks < p_value_cutoff)[0]
    significant_diff_reactions = [cobra_reactions_str[i] for i in significant_diff_indices]
    
    not_significant_diff_indices = np.where(corrected_p_values_ks >= p_value_cutoff)[0]
    not_significant_diff_reactions = [cobra_reactions_str[i] for i in not_significant_diff_indices]
    
    return significant_diff_indices, significant_diff_reactions, not_significant_diff_reactions