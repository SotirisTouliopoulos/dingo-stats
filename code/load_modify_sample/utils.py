
from cobra.io import read_sbml_model
from dingo import MetabolicNetwork, PolytopeSampler
from cobra.sampling.optgp import OptGPSampler
from cobra.util.solver import linear_reaction_coefficients
from dingo import set_default_solver
from gapsplit import gapsplit
from cobra.flux_analysis.loopless import add_loopless


def load_model(filepath):
    cobra_model = read_sbml_model(filepath)
    cobra_reactions = cobra_model.reactions
    
    dingo_model = MetabolicNetwork.from_cobra_model(cobra_model)
    dingo_reactions = dingo_model.reactions

    return cobra_model, cobra_reactions, dingo_model, dingo_reactions


def get_objective_functions(cobra_model):
    objectives_dict = linear_reaction_coefficients(cobra_model)
    objective_functions = list(objectives_dict.keys())
    objective_functions_ids = [rxn.id for rxn in objective_functions]
    
    return objective_functions_ids


def get_reaction_bounds(cobra_model):
    reaction_bounds_dict = { }
    reaction_ids = [ reaction.id for reaction in cobra_model.reactions ]
                
    for reaction_id in reaction_ids:
        bounds = cobra_model.reactions.get_by_id(reaction_id).bounds
        lb = round(bounds[0], 6)
        up = round(bounds[1], 6)
        
        reaction_bounds_dict[reaction_id] = (lb, up)
        
    return reaction_bounds_dict


def modify_model(cobra_model, objective_function="", optimal_percentage=100, reaction_bounds={}):
    modified_cobra_model = cobra_model.copy()
    
    # first define reaction bounds from dictionary
    if len(reaction_bounds) >= 1:
        for reaction, rbounds in reaction_bounds.items():
            rxn = modified_cobra_model.reactions.get_by_id(reaction)
            rxn.lower_bound, rxn.upper_bound = rbounds
        
    # then modify object function and optimal percentage (this may overwrite bounds for objective function from previous dictionary)
    modified_cobra_model.objective = objective_function
    modified_cobra_model.reactions.get_by_id(objective_function).upper_bound = 1000
    
    fba_solution = modified_cobra_model.optimize()
    objective = fba_solution.objective_value
    
    fraction_of_optimum = float(optimal_percentage / 100)
    # in case due to numeric errors it gets a higher value
    if fraction_of_optimum > 1:
        fraction_of_optimum = 1
    
    # in case due to numeric errors it gets a higher value
    minimum_objective_value = round( (objective * fraction_of_optimum), 6)
    if minimum_objective_value > objective:
        minimum_objective_value = minimum_objective_value
        
    # to avoid errors
    if optimal_percentage == 100:
        minimum_objective_value = minimum_objective_value - 1e-3
    
    modified_cobra_model.reactions.get_by_id(objective_function).lower_bound = minimum_objective_value
    
    modified_dingo_model = MetabolicNetwork.from_cobra_model(modified_cobra_model)
    modified_dingo_model.set_opt_percentage(optimal_percentage)
    
    return modified_cobra_model, modified_dingo_model
    

def sample_optgp(cobra_model, n_samples = 3000, reaction_in_rows = True):
    sampler_optgp = OptGPSampler(cobra_model)
    samples_optgp = sampler_optgp.sample(n=n_samples)
    samples_optgp = samples_optgp.to_numpy()
    
    if reaction_in_rows == True:
        samples_optgp = samples_optgp.T
        
    return samples_optgp
    
    
def sample_dingo(dingo_model, reaction_in_rows = True, ess=1000, psrf = False):
    set_default_solver("gurobi")

    sampler_dingo_default = PolytopeSampler(dingo_model)
    samples_dingo_default = sampler_dingo_default.generate_steady_states(ess=ess, psrf=psrf)
    
    if reaction_in_rows == True:
        pass
    else:
        samples_dingo_default = samples_dingo_default.T
        
    return samples_dingo_default


def sample_gapsplit(cobra_model, n_samples = 3000, reaction_in_rows = True, add_loopless_cobrapy = False):
    
    ec_cobra_model = cobra_model.copy()
    
    if add_loopless_cobrapy == True:
        add_loopless(ec_cobra_model)
    
    gapsplit_samples = gapsplit(
        ec_cobra_model, 
        n=n_samples, 
        gurobi_direct=False,
        fraction_of_optimum=0.5
    )
    
    gapsplit_samples = gapsplit_samples.to_numpy()
    
    if reaction_in_rows == True:
        gapsplit_samples = gapsplit_samples.T
        
    return gapsplit_samples