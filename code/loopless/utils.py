
from cobra.io import read_sbml_model
from dingo import MetabolicNetwork, PolytopeSampler
from cobra.flux_analysis import flux_variability_analysis
from cobra.sampling.optgp import OptGPSampler
import cobra
from cobra.util.solver import linear_reaction_coefficients
from dingo import set_default_solver
import numpy as np
from scipy import stats
import pickle
from cobra.flux_analysis.loopless import loopless_solution
import pandas as pd
import matplotlib.pyplot as plt


def load_model(filepath):
    _cobra_model = read_sbml_model(filepath)
    _cobra_reactions = _cobra_model.reactions
    
    _dingo_model = MetabolicNetwork.from_cobra_model(_cobra_model)
    _dingo_reactions = _dingo_model.reactions

    return _cobra_model, _cobra_reactions, _dingo_model, _dingo_reactions


def get_objective_functions(cobra_model):
    _objectives_dict = linear_reaction_coefficients(cobra_model)
    _objective_functions = list(_objectives_dict.keys())
    _objective_functions_ids = [rxn.id for rxn in _objective_functions]
    
    return _objective_functions_ids


def get_reaction_bounds(cobra_model):
    _reaction_bounds_dict = { }
    _reaction_ids = [ reaction.id for reaction in cobra_model.reactions ]
                
    for reaction_id in _reaction_ids:
        _bounds = cobra_model.reactions.get_by_id(reaction_id).bounds
        _lb = round(_bounds[0], 6)
        _up = round(_bounds[1], 6)
        
        _reaction_bounds_dict[reaction_id] = (_lb, _up)
        
    return _reaction_bounds_dict


def modify_model(cobra_model, objective_function="", optimal_percentage=100, reaction_bounds={}):
    _modified_cobra_model = cobra_model.copy()
    
    # first define reaction bounds from dictionary
    if len(reaction_bounds) >= 1:
        for _reaction, _rbounds in reaction_bounds.items():
            _rxn = _modified_cobra_model.reactions.get_by_id(_reaction)
            _rxn.lower_bound, _rxn.upper_bound = _rbounds
        
    # then modify object function and optimal percentage (this may overwrite bounds for objective function from previous dictionary)
    _modified_cobra_model.objective = objective_function
    _modified_cobra_model.reactions.get_by_id(objective_function).upper_bound = 1000
    
    _fba_solution = _modified_cobra_model.optimize()
    _objective = _fba_solution.objective_value
    
    _fraction_of_optimum = float(optimal_percentage / 100)
    # in case due to numeric errors it gets a higher value
    if _fraction_of_optimum > 1:
        _fraction_of_optimum = 1
    
    # in case due to numeric errors it gets a higher value
    _minimum_objective_value = round( (_objective * _fraction_of_optimum), 6)
    if _minimum_objective_value > _objective:
        _minimum_objective_value = _minimum_objective_value
        
    # to avoid errors
    if optimal_percentage == 100:
        _minimum_objective_value = _minimum_objective_value - 1e-3
    
    _modified_cobra_model.reactions.get_by_id(objective_function).lower_bound = _minimum_objective_value
    
    _modified_dingo_model = MetabolicNetwork.from_cobra_model(_modified_cobra_model)
    _modified_dingo_model.set_opt_percentage(optimal_percentage)
    
    return _modified_cobra_model, _modified_dingo_model
    

def fva_dictionary(cobra_model, optimal_percentage = 100, loopless = False):
    _fraction_of_optimum = float(optimal_percentage / 100)

    _fva = flux_variability_analysis(cobra_model, reaction_list=cobra_model.reactions, fraction_of_optimum=_fraction_of_optimum, loopless=loopless)
    _fva_dict = _fva.apply(lambda row: (row['minimum'], row['maximum']), axis=1).to_dict()
    
    return _fva_dict


def sample_optgp(cobra_model, n_samples = 3000, reaction_in_rows = True):
    _sampler_optgp = OptGPSampler(cobra_model)
    _samples_optgp = _sampler_optgp.sample(n=n_samples)
    _samples_optgp = _samples_optgp.to_numpy()
    
    if reaction_in_rows == True:
        _samples_optgp = _samples_optgp.T
        
    return _samples_optgp
    
    
def sample_dingo(dingo_model, reaction_in_rows = True, ess=1000, psrf = False):
    set_default_solver("gurobi")

    _sampler_dingo_default = PolytopeSampler(dingo_model)
    _samples_dingo_default = _sampler_dingo_default.generate_steady_states(ess=ess, psrf=psrf)
    
    if reaction_in_rows == True:
        pass
    else:
        _samples_dingo_default = _samples_dingo_default.T
        
    return _samples_dingo_default


def export_to_pickle(samples, filename):
    with open(filename, "wb") as _hopsy_samples_file: 
        pickle.dump(samples, _hopsy_samples_file)


def sampling_statistics(samples, reactions_ids_list=None, reaction_id=""):
    _reaction_index = reactions_ids_list.index("FRD7")
    
    _mean = np.mean(samples[_reaction_index])
    _min = np.min(samples[_reaction_index]) 
    _max = np.max(samples[_reaction_index])
    _std = np.std(samples[_reaction_index])
    _skewness = stats.skew(samples[_reaction_index])
    _kurtosis = stats.kurtosis(samples[_reaction_index])
    
    return _mean, _min, _max, _std, _skewness, _kurtosis
        
    
def get_loopless_solutions_from_samples(samples, cobra_model):
    _cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]
    _samples = samples.copy()
        
    # if provided sampling dataset has reactions as rows ==> transpose
    if _samples.shape[0] == len(_cobra_reactions_str):
        _samples = _samples.T 
    
    _loopless_solutions_pandas_series_ec_default = []
    for i in range((_samples).shape[0]):
        
        _sample = (_samples)[i]
        _sample_reactions_dictionary = {_k:_v for _k,_v in zip(_cobra_reactions_str, _sample)}
        _loopless_sample = loopless_solution(model=cobra_model, fluxes=_sample_reactions_dictionary)
        
        _loopless_solutions_pandas_series_ec_default.append(_loopless_sample.fluxes)

    df = pd.concat(_loopless_solutions_pandas_series_ec_default, axis=1)
    _samples_ec_model_loopless_solutions = df.to_numpy()

    return _samples_ec_model_loopless_solutions


def calculate_affected_samples(initial_samples, loopless_samples, cobra_model, tol_reaction_difference=0.01, tol_reactions_count=1):
    _cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]
    _initial_samples = initial_samples.copy()
    _loopless_samples = loopless_samples.copy()
        
    # if provided sampling dataset has reactions as rows ==> transpose
    if _initial_samples.shape[0] == len(_cobra_reactions_str):
        _initial_samples = _initial_samples.T
        
    # if provided sampling dataset has reactions as rows ==> transpose
    if _loopless_samples.shape[0] == len(_cobra_reactions_str):
        _loopless_samples = _loopless_samples.T
        
    # find the ranges of each reaction to normalize the reaction difference tolerance
    _reactions_ranges = {_cobra_reactions_str[i]: _initial_samples[:, i].max() - _initial_samples[:, i].min() for i in range(_initial_samples.shape[1])}


    _affected_reactions_count = []
    _total_affected_samples = 0
    for i in range((_initial_samples).shape[0]):
        
        df = pd.DataFrame({
        'value_1': _initial_samples[i],
        'value_2': _loopless_samples[i],
        'difference': abs(_initial_samples[i] - _loopless_samples[i])
        })
        df['norm_difference'] = df['difference'] / list(_reactions_ranges.values())


        df_filtered = df[df['norm_difference'] > tol_reaction_difference]
        _affected_reactions_count.append(df_filtered.shape[0])
        
        if df_filtered.shape[0] >= tol_reactions_count:
            _total_affected_samples += 1
    
    return _affected_reactions_count, _total_affected_samples


def plot_grid_95_reactions(samples, cobra_model, nrows=20, ncols=5):
    _cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]
    _samples = samples.copy()
    
    if nrows * ncols < len(_cobra_reactions_str):
        raise Exception("Change dimensions of the plot (nrows, ncols) to store all distributions")
    
    # if provided sampling dataset has reactions as rows ==> transpose
    if _samples.shape[0] == len(_cobra_reactions_str):
        _samples = _samples.T

    _df = pd.DataFrame(_samples)
    _df.columns = _cobra_reactions_str

    _nrows = nrows
    _ncols = ncols
    _fig, _axes = plt.subplots(nrows=_nrows, ncols=_ncols, figsize=(4*_ncols, 3*_nrows))
    _axes = np.array(_axes).flatten()

    # Plot each column in its subplot
    for i, col in enumerate(_df.columns):
        ax = _axes[i]
        ax.hist(_df[col], bins=20, color='skyblue', edgecolor='black')
        ax.set_title(col)
        ax.grid(True)

    # If there are unused axes (e.g., more grid spots than columns), hide them
    for j in range(len(_df.columns), len(_axes)):
        _fig.delaxes(_axes[j])  # or: axes[j].axis('off')

    plt.tight_layout()
    plt.show()
    
    
def set_bounds_from_loopless_solution_samples(loopless_solution_samples, cobra_model):
    _modified_cobra_model = cobra_model.copy()
    _cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]
    _initial_samples = loopless_solution_samples.copy()
        
    # if provided sampling dataset has reactions as rows ==> transpose
    if _initial_samples.shape[0] == len(_cobra_reactions_str):
        _initial_samples = _initial_samples.T
        
    # find the ranges of each reaction to normalize the reaction difference tolerance
    _reaction_bounds = {
        _cobra_reactions_str[i]: (np.min(_initial_samples[:, i]), np.max(_initial_samples[:, i]))
        for i in range(_initial_samples.shape[1])
    }

    for _reaction, _rbounds in _reaction_bounds.items():
        
        _rxn = _modified_cobra_model.reactions.get_by_id(_reaction)
        _rxn_lb = round(_rbounds[0], 6)
        _rxn_up = round(_rbounds[1], 6)

        _rxn.lower_bound = _rxn_lb
        _rxn.upper_bound = _rxn_up
            
    return _modified_cobra_model
