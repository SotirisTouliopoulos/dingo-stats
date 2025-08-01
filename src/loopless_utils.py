
import numpy as np
from cobra.flux_analysis.loopless import loopless_solution
import pandas as pd
import plotly.express as px
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis import pfba
from typing import Dict, Tuple, List
from cobra import Model, Reaction
from numpy.typing import NDArray

def loops_enumeration_from_fva(
    ec_cobra_model: Model, 
    fraction_of_optimum: float = 1.0
) -> List:
    """
    Function that aims to identify loopy reaction from Flux Variability Analysis.
    This method was introduced here: https://doi.org/10.1016/j.bpj.2010.12.3707
    
    Keyword arguments:
    cobra_model (Model) -- cobra model object
    fraction_of_optimum (float) -- Float variable that defines the `fraction_of_optimum` parameter used in `flux_variability_analysis`

    Returns:
    loopy_reactions (List) -- List with reactions classified as loopy based on the algorithm
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
    
    loopy_reactions = list(zip(significant_diff.index.tolist(), significant_diff["span_difference"].tolist()))

    return loopy_reactions


def loopy_reactions_turned_off_in_pfba(
    loopy_cobra_model: Model, 
    loopy_reactions: List, 
    tol: float = 1e-6
) -> List:
    """
    Function that identifies which reactions from the loopy reactions set calculated with the `loops_enumeration_from_fva` function 
    can be turned off when performing a `pFBA` or a `loopless_solution` to `pFBA` results

    Keyword arguments:
    loopy_cobra_model (Model) -- cobra model object possibly containing loops (default model without any preproccess)
    loopy_reactions (List) -- loopy_reactions (List) -- List with reactions classified as loopy
    tol (float) -- Tolerance value used to classify a numeric change as important or not

    Returns:
    turned_off_reactions (List) -- List with reactions turned of in `pFBA` with `loopless_solution` applied 
    """

    fba_model = loopy_cobra_model.copy()
    solution_pfba = pfba(fba_model)

    pfba_loopless_solution = loopless_solution(loopy_cobra_model, solution_pfba)

    turned_off_reactions = []
    for rxn in loopy_reactions:
        flux_pfba_loopless = pfba_loopless_solution[rxn]
        
        if abs(flux_pfba_loopless) < tol:
            turned_off_reactions.append(rxn)
            
    return turned_off_reactions


def get_loopless_solutions_from_samples(
    samples: NDArray[np.float64],
    cobra_model: Model
) -> NDArray[np.float64]:
    """
    Function that calculates the `loopless_solution` from each sample and saves results into a new numpy 2D array
    
    Keyword arguments:
    samples (NDArray[np.float64]) -- Numpy 2D array of the samples
    cobra_model (model) -- cobra model object
    
    Returns:
    samples_loopless_solutions (NDArray[np.float64]) -- Numpy 2D array of the samples after application of `loopless_solution`
    """
    
    cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]
    samples = samples.copy()
    
    if samples.shape[0] == samples.shape[1]:
        raise ValueError("Samples array provided has equal rows and columns dimensions. Please change the number of samples")
        
    # if provided sampling dataset has reactions as rows ==> transpose
    if samples.shape[0] == len(cobra_reactions_str):
        samples = samples.T 
    
    loopless_solutions_pandas_series_default = []
    for i in range((samples).shape[0]):
        
        sample = (samples)[i]
        sample_reactions_dictionary = {k:v for k,v in zip(cobra_reactions_str, sample)}
        loopless_sample = loopless_solution(model=cobra_model, fluxes=sample_reactions_dictionary)
        
        loopless_solutions_pandas_series_default.append(loopless_sample.fluxes)

    df = pd.concat(loopless_solutions_pandas_series_default, axis=1)
    samples_loopless_solutions = df.to_numpy()

    return samples_loopless_solutions


def calculate_affected_samples(
    initial_samples: NDArray[np.float64], 
    loopless_samples: NDArray[np.float64],
    cobra_model: Model,
    tol_reaction_difference: float = 0.01,
    tol_reactions_count: int = 1
) -> Tuple[List, int]:
    """
    Function that given the default and the `get_loopless_solutions_from_samples` samples identifies how many samples have significant differences.
    This is done by calculating (for each sample) the distance (change) of the reactions from the given arrays.
    Reactions that have a difference higher than a given threshold are counted and the corresponding sample can potentially be classified as a loopy sample.
    However, if user does not consider importance when only 1 reaction overcomes this tolerance, he can specify the least number 
    of reactions needed for this classification (e.g 2, 5, 10 ..).
    
    Keyword arguments:
    initial_samples (NDArray[np.float64]) -- Numpy 2D array of the default samples
    loopless_samples (NDArray[np.float64]) -- Numpy 2D array of the default samples after application of the `get_loopless_solutions_from_samples`
    cobra_model (Model) -- ocbra model object
    tol_reaction_difference (float) -- Tolerance value to classify the difference in reactions before and after `get_loopless_solutions_from_samples` as important
    tol_reactions_count (int) -- Tolerance value to classify a sample as loopy based on the least amount of significantly altered reactions
    
    Returns:
    Tuple[List, int]
        affected_reactions_count (List) -- List containing the number of affected (significantly altered) reactions across all samples
        total_affected_samples (int) -- The total number of affected (significantly altered) samples
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


def calculate_distances_from_samples(
    initial_samples: NDArray[np.float64], 
    loopless_samples: NDArray[np.float64], 
    cobra_model: Model
) -> NDArray[np.float64]:
    """
    Function that calculates the euclidean distance between a sampling dataset before and after applying the `get_loopless_solutions_from_samples`.
    Distances are calculated between different samples, so if user has 3000 samples, he will end up with 3000 distances.
    
    Keyword arguments:
    initial_samples (NDArray[np.float64]) -- Numpy 2D array of the default samples
    loopless_samples (NDArray[np.float64]) -- Numpy 2D array of the default samples after application of the `get_loopless_solutions_from_samples`
    cobra_model (Model) -- cobra model object
    
    Returns:
    distances_array (NDArray[np.float64]) -- Numpy 1D array of the samples' euclidean distances before and after the `get_loopless_solutions_from_samples`
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


def calculate_distances_from_reactions(
    initial_samples: NDArray[np.float64], 
    loopless_samples: NDArray[np.float64], 
    cobra_model: Model
) -> NDArray[np.float64]:
    """
    Function that calculates the euclidean distance between a sampling dataset before and after applying the `get_loopless_solutions_from_samples`.
    Distances are calculated between different reactions, so if user provides samples of 100 reactions, he will end up with 100 distances.    
    
    Keyword arguments:
    initial_samples (NDArray[np.float64]) -- Numpy 2D array of the default samples
    loopless_samples (NDArray[np.float64]) -- Numpy 2D array of the default samples after application of the `get_loopless_solutions_from_samples`
    cobra_model (Model) -- cobra model object
    
    Returns:
    distances_array (NDArray[np.float64]) -- Numpy 1D array of the reactions' euclidean distances before and after the `get_loopless_solutions_from_samples`
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


def violin_plot_samples_distances(
    distances_array: NDArray[np.float64], 
    title: str = "", 
    exponentformat: str = "none"):
    """
    Function that creates a violin plot from a given `distances_array`
    
    Keyword arguments:
        distances_array (NDArray[np.float64]) -- Numpy 1D array of either the reactions' or samples' euclidean distances before and after 
                                                 the `get_loopless_solutions_from_samples`
        title (str) -- Title for the plot
        exponentformat (str) -- Parameter used to determine the numeric scientific notation in the plot (determines a formatting rule for the tick exponents). 
                                Available options: "e", "E", "power", "SI", "B"
    """
    
    df = pd.DataFrame({'Distance': distances_array})
    fig = px.violin(df, y='Distance', box=True, points='all', title=title)
    
    fig.update_layout(
    yaxis=dict(
        exponentformat=exponentformat)
    )
    
    fig.show()