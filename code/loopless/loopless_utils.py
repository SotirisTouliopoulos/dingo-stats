
import numpy as np
from cobra.flux_analysis.loopless import loopless_solution
import pandas as pd
        
    
def get_loopless_solutions_from_samples(samples, cobra_model):
    cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]
    samples = samples.copy()
        
    # if provided sampling dataset has reactions as rows ==> transpose
    if samples.shape[0] == len(cobra_reactions_str):
        samples = samples.T 
    
    loopless_solutions_pandas_series_ec_default = []
    for i in range((samples).shape[0]):
        
        sample = (samples)[i]
        sample_reactions_dictionary = {k:v for k,v in zip(cobra_reactions_str, sample)}
        loopless_sample = loopless_solution(model=cobra_model, fluxes=sample_reactions_dictionary)
        
        loopless_solutions_pandas_series_ec_default.append(loopless_sample.fluxes)

    df = pd.concat(loopless_solutions_pandas_series_ec_default, axis=1)
    samples_ec_model_loopless_solutions = df.to_numpy()

    return samples_ec_model_loopless_solutions


def calculate_affected_samples(initial_samples, loopless_samples, cobra_model, tol_reaction_difference=0.01, tol_reactions_count=1):
    cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]
    initial_samples = initial_samples.copy()
    loopless_samples = loopless_samples.copy()
        
    # if provided sampling dataset has reactions as rows ==> transpose
    if initial_samples.shape[0] == len(cobra_reactions_str):
        initial_samples = initial_samples.T
        
    # if provided sampling dataset has reactions as rows ==> transpose
    if loopless_samples.shape[0] == len(cobra_reactions_str):
        loopless_samples = loopless_samples.T
        
    # find the ranges of each reaction to normalize the reaction difference tolerance
    reactions_ranges = {cobra_reactions_str[i]: initial_samples[:, i].max() - initial_samples[:, i].min() for i in range(initial_samples.shape[1])}


    affected_reactions_count = []
    total_affected_samples = 0
    for i in range((initial_samples).shape[0]):
        
        df = pd.DataFrame({
        'value_1': initial_samples[i],
        'value_2': loopless_samples[i],
        'difference': abs(initial_samples[i] - loopless_samples[i])
        })
        df['norm_difference'] = df['difference'] / list(reactions_ranges.values())


        df_filtered = df[df['norm_difference'] > tol_reaction_difference]
        affected_reactions_count.append(df_filtered.shape[0])
        
        if df_filtered.shape[0] >= tol_reactions_count:
            total_affected_samples += 1
    
    return affected_reactions_count, total_affected_samples
    
    
def set_bounds_from_loopless_solution_samples(loopless_solution_samples, cobra_model):
    modified_cobra_model = cobra_model.copy()
    cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]
    initial_samples = loopless_solution_samples.copy()
        
    # if provided sampling dataset has reactions as rows ==> transpose
    if initial_samples.shape[0] == len(cobra_reactions_str):
        initial_samples = initial_samples.T
        
    # find the ranges of each reaction to normalize the reaction difference tolerance
    reaction_bounds = {
        cobra_reactions_str[i]: (np.min(initial_samples[:, i]), np.max(initial_samples[:, i]))
        for i in range(initial_samples.shape[1])
    }

    for reaction, rbounds in reaction_bounds.items():
        
        rxn = modified_cobra_model.reactions.get_by_id(reaction)
        rxn_lb = round(rbounds[0], 6)
        rxn_up = round(rbounds[1], 6)

        rxn.lower_bound = rxn_lb
        rxn.upper_bound = rxn_up
            
    return modified_cobra_model


def calculate_distances_from_samples(initial_samples, loopless_samples, cobra_model):
    _cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]
    _initial_samples = initial_samples.copy()
    _loopless_samples = loopless_samples.copy()
        
    # if provided sampling dataset has reactions as rows ==> transpose
    if _initial_samples.shape[0] == len(_cobra_reactions_str):
        _initial_samples = _initial_samples.T
        
    # if provided sampling dataset has reactions as rows ==> transpose
    if _loopless_samples.shape[0] == len(_cobra_reactions_str):
        _loopless_samples = _loopless_samples.T
        
    _conditions = [ _initial_samples, _loopless_samples ]
    _selected_comparisons = [(0, 1)]

    distances_array = np.zeros(len(_initial_samples))
    for row in range(len(_initial_samples)):
        for i, j in _selected_comparisons:
            
            distances_array[row] = np.linalg.norm(_conditions[i][row] - _conditions[j][row])

    return distances_array


def calculate_distances_from_reactions(initial_samples, loopless_samples, cobra_model):
    _cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]
    _initial_samples = initial_samples.copy()
    _loopless_samples = loopless_samples.copy()
        
    # if provided sampling dataset has reactions as cols ==> transpose
    if _initial_samples.shape[1] == len(_cobra_reactions_str):
        _initial_samples = _initial_samples.T
        
    # if provided sampling dataset has reactions as cols ==> transpose
    if _loopless_samples.shape[1] == len(_cobra_reactions_str):
        _loopless_samples = _loopless_samples.T
        
    _conditions = [ _initial_samples, _loopless_samples ]
    _selected_comparisons = [(0, 1)]

    distances_array = np.zeros(len(_cobra_reactions_str))
    for row in range(len(_cobra_reactions_str)):
        for i, j in _selected_comparisons:
            
            distances_array[row] = np.linalg.norm(_conditions[i][row] - _conditions[j][row])

    return distances_array
