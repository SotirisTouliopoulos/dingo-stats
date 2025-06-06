

from statsmodels.stats.multitest import multipletests
from scipy import stats
import numpy as np
        

def significantly_altered_reactions(conditions=[], selected_comparisons = [(0, 1)], 
                                    cobra_model=None, p_value_cutoff=0.05, fold_change_cutoff=1.2, std_cutoff=1e-2):
    
    """
    Function that takes as input at least 2 conditions to compare (diffeent sampling dataset) 
    and performs the KS non-parametric test and based on a multiple-comparisons corrected p-value and fold change classifies
    reactions as significantly altered or not. This function can provide reactions that can further filter graphs
    """
    
    cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]
        
    # if provided sampling dataset has reactions as cols ==> transpose
    for condition in conditions:
        if condition.shape[1] == len(cobra_reactions_str):
            condition = condition.T
    
    # Store p-values for each row for KS test
    p_values = {row: [] for row in range(len(cobra_reactions_str)) }
    # Store ks-values for each row for KS test
    ks_values = {row: [] for row in range(len(cobra_reactions_str)) }
    # Store fold-change-values for each row for KS test
    fold_change_values = {row: [] for row in range(len(cobra_reactions_str)) }

    for row in range(len(cobra_reactions_str)):
        for i, j in selected_comparisons:
            
            if (np.std(conditions[i][row]) > std_cutoff) or (np.std(conditions[j][row]) > std_cutoff):
            
                ks, p = stats.ks_2samp(conditions[i][row], conditions[j][row], alternative='two-sided')
                p_values[row].append(p)
                ks_values[row].append(ks)
                
                fold_change = np.absolute( (np.mean(conditions[i][row]) - np.mean(conditions[j][row]) ) / (np.mean(conditions[j][row]) + 1e-8) )
                fold_change_values[row].append(fold_change)
                
            else:
                p_values[row].append(1)
                ks_values[row].append(1)
                fold_change_values[row].append(0)


    p_values_copy = list(p_values.values())
    flat_p_values = np.array(p_values_copy).flatten()
    
    fold_change_values = np.array(list(fold_change_values.values()))


    # Apply FDR correction
    _, corrected_p_values, _, _ = multipletests(flat_p_values, method='fdr_bh')

    # Reshape corrected p-values back to the original matrix shape
    corrected_p_values = corrected_p_values.reshape(np.array(p_values_copy).shape)
    
    significant_diff_indices = np.where(
    np.logical_and(
        corrected_p_values < p_value_cutoff,
        fold_change_values > fold_change_cutoff
    ) )[0]
    
    significant_diff_reactions = [cobra_reactions_str[i] for i in significant_diff_indices]
    
    not_significant_diff_indices = np.where(corrected_p_values >= p_value_cutoff)[0]
    not_significant_diff_reactions = [cobra_reactions_str[i] for i in not_significant_diff_indices]
    
    return significant_diff_indices, significant_diff_reactions, not_significant_diff_reactions