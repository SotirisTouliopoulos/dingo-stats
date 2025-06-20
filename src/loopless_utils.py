
import numpy as np
from cobra.flux_analysis.loopless import loopless_solution
import pandas as pd
import plotly.express as px
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis import pfba


def loops_enumeration_from_fva(ec_cobra_model, fraction_of_optimum=1):
    """
    Function that aims to identify loopy reaction from FVA.
    This method was introduced here: https://doi.org/10.1016/j.bpj.2010.12.3707
    """
    
    fva_noloop = flux_variability_analysis(ec_cobra_model, loopless=False, fraction_of_optimum=fraction_of_optimum)

    # Perform FVA with loopless constraints
    fva_loopless = flux_variability_analysis(ec_cobra_model, loopless=True, fraction_of_optimum=fraction_of_optimum)

    # Calculate span = max - min for each
    fva_noloop['span_noloop'] = fva_noloop['maximum'] - fva_noloop['minimum']
    fva_loopless['span_loopless'] = fva_loopless['maximum'] - fva_loopless['minimum']

    # Combine the results into a single DataFrame
    fva_comparison = pd.DataFrame({
        'span_noloop': fva_noloop['span_noloop'],
        'span_loopless': fva_loopless['span_loopless']
    })

    # Calculate the difference in spans (optional)
    fva_comparison['span_difference'] = fva_comparison['span_noloop'] - fva_comparison['span_loopless']
    
    significant_diff = fva_comparison[fva_comparison['span_difference'] > 1e-6]
    
    significant_diff = list(zip(significant_diff.index.tolist(), significant_diff["span_difference"].tolist()))

    return significant_diff


def loopy_reactions_turned_off_in_pfba(loopy_cobra_model, loopy_reactions):
    """
    Identify which subset of loopy reactions calculated with the above function 
    can be turned off when performing a pFBA or a loopless solutions to pFBA results
    """

    fba_model = loopy_cobra_model.copy()
    pfba_solution = pfba(fba_model)

    lloopless_solution = loopless_solution(loopy_cobra_model)

    lloopless_pfba_sol = loopless_solution(loopy_cobra_model, pfba_solution)

    turned_off_in_pfba = []
    tol = 1e-6

    for rxn in loopy_reactions:
        
        flux_pfba = pfba_solution.fluxes[rxn]
        flux_loopless = lloopless_solution.fluxes[rxn]
        flux_pfba_loopless = lloopless_pfba_sol[rxn]

        if abs(flux_pfba_loopless) < tol:

            turned_off_in_pfba.append(rxn)
            
    return turned_off_in_pfba


def get_loopless_solutions_from_samples(samples, cobra_model):
    """
    Function that calculate the loopless solution from each sample and saves results into a new numpy 2D array
    """
    
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
    """
    Function that given the default and the lopless solutions samples identifies how many samples have significant differences.
    This is done by calculating (for each sample) the distance across the reactions between the default and loopless solution.
    Reactions that have a difference higher than a given threshold are counted and sample can potentially be classified as loopy.
    However, if user does not assume a sample as loopy when only 1 reaction overcomes this value, he can specify the least number 
    of reactions needed for this classification. 
    """
    
    cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]
    initial_samples = initial_samples.copy()
    loopless_samples = loopless_samples.copy()
        
    # if provided sampling dataset has reactions as rows ==> transpose
    if initial_samples.shape[0] == len(cobra_reactions_str):
        initial_samples = initial_samples.T
        
    # if provided sampling dataset has reactions as rows ==> transpose
    if loopless_samples.shape[0] == len(cobra_reactions_str):
        loopless_samples = loopless_samples.T
        
    # find the flux range of each reaction from the samples in order to normalize the reaction difference tolerance
    reactions_ranges = {cobra_reactions_str[i]: initial_samples[:, i].max() - initial_samples[:, i].min() for i in range(initial_samples.shape[1])}


    affected_reactions_count = []
    total_affected_samples = 0
    # loop through every sample and find significantly altered reactions
    for i in range((initial_samples).shape[0]):
        
        df = pd.DataFrame({
        'value_1': initial_samples[i],
        'value_2': loopless_samples[i],
        'difference': abs(initial_samples[i] - loopless_samples[i])
        })
        df['norm_difference'] = df['difference'] / list(reactions_ranges.values())

        # keep only reactions over the difference threshold
        df_filtered = df[df['norm_difference'] > tol_reaction_difference]
        affected_reactions_count.append(df_filtered.shape[0])
        
        # user can classify a sample as affected when at least X reactions are removed (default 1 reaction)
        if df_filtered.shape[0] >= tol_reactions_count:
            total_affected_samples += 1
    
    return affected_reactions_count, total_affected_samples


def calculate_distances_from_samples(initial_samples, loopless_samples, cobra_model):
    """
    Function that calculates the euclidean distance between a sampling dataset and itself after applyig the loopless solutions.
    Distances are calculated between different samples, so if you have 3000 samples you will end up with 3000 distances.
    """
    
    cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]
    initial_samples = initial_samples.copy()
    loopless_samples = loopless_samples.copy()
        
    # if provided sampling dataset has reactions as rows ==> transpose
    if initial_samples.shape[0] == len(cobra_reactions_str):
        initial_samples = initial_samples.T
        
    # if provided sampling dataset has reactions as rows ==> transpose
    if loopless_samples.shape[0] == len(cobra_reactions_str):
        loopless_samples = loopless_samples.T
        
    conditions = [ initial_samples, loopless_samples ]
    selected_comparisons = [(0, 1)]

    distances_array = np.zeros(len(initial_samples))
    for row in range(len(initial_samples)):
        for i, j in selected_comparisons:
            
            distances_array[row] = np.linalg.norm(conditions[i][row] - conditions[j][row])

    return distances_array


def calculate_distances_from_reactions(initial_samples, loopless_samples, cobra_model):
    """
    Same as above but distances are calculated between reactions.
    """
    
    cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]
    initial_samples = initial_samples.copy()
    loopless_samples = loopless_samples.copy()
        
    # if provided sampling dataset has reactions as cols ==> transpose
    if initial_samples.shape[1] == len(cobra_reactions_str):
        initial_samples = initial_samples.T
        
    # if provided sampling dataset has reactions as cols ==> transpose
    if loopless_samples.shape[1] == len(cobra_reactions_str):
        loopless_samples = loopless_samples.T
        
    conditions = [ initial_samples, loopless_samples ]
    selected_comparisons = [(0, 1)]

    distances_array = np.zeros(len(cobra_reactions_str))
    for row in range(len(cobra_reactions_str)):
        for i, j in selected_comparisons:
            
            distances_array[row] = np.linalg.norm(conditions[i][row] - conditions[j][row])

    return distances_array


def violin_plot_samples_distances(distances_array, title=""):
    """
    Function that creates a violin plot from a distances_array (this array is calculated from the above 2 functions)
    """
    
    df = pd.DataFrame({'Distance': distances_array})
    fig = px.violin(df, y='Distance', box=True, points='all', title=title)
    fig.show()