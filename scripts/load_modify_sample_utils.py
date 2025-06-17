
from cobra.io import read_sbml_model
from dingo import MetabolicNetwork, PolytopeSampler
from cobra.sampling.optgp import OptGPSampler
from cobra.util.solver import linear_reaction_coefficients
from dingo import set_default_solver
from gapsplit import gapsplit
from cobra.flux_analysis.loopless import add_loopless
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import pickle
import logging


def load_model(filepath):
    """
    Function to load and return a model in both cobra and dingo formats from given path.
    It additionally returns 2 lists, one with the corresponding model reactions (cobra_reactions) 
    and one with their IDs (dingo_reactions)
    """

    cobra_model = read_sbml_model(filepath)
    cobra_reactions = cobra_model.reactions

    dingo_model = MetabolicNetwork.from_cobra_model(cobra_model)
    dingo_reactions = dingo_model.reactions

    return cobra_model, cobra_reactions, dingo_model, dingo_reactions


def get_objective_functions(cobra_model):
    """
    Function to find and return the objective functions of a cobra model
    """

    objectives_dict         = linear_reaction_coefficients(cobra_model)
    objective_functions     = list(objectives_dict.keys())
    objective_functions_ids = [rxn.id for rxn in objective_functions]

    return objective_functions_ids


def get_reaction_bounds(cobra_model):
    """
    Function to create and return a dictionary of a cobra's model reaction IDs as keys and 
    their corresponding bounds as values
    """

    reaction_bounds_dict = {}
    reaction_ids = [reaction.id for reaction in cobra_model.reactions]

    for reaction_id in reaction_ids:
        bounds = cobra_model.reactions.get_by_id(reaction_id).bounds
        lb = round(bounds[0], 6)
        up = round(bounds[1], 6)

        reaction_bounds_dict[reaction_id] = (lb, up)

    return reaction_bounds_dict


def modify_model(
    cobra_model, objective_function=None, optimal_percentage=100, reaction_bounds={}, objective_direction="max"
):
    """
    Function that modifies a given cobra model.

    Users can change objective function, optimal percentage, define custom reaction bounds and change objective 
    direction.

    NOTE (Haris Zafeiropoulos, 2025-06-16):
    Key-feature, lower bound of the objective will be fixed based on the optimal percentage, in all cases.. 
    """

    modified_cobra_model = cobra_model.copy()

    # Reaction bounds from a user-provided dictionary,
    # where keys are reaction ids and values the min/max flux bounds
    if len(reaction_bounds) >= 1:
        for reaction, rbounds in reaction_bounds.items():
            rxn = modified_cobra_model.reactions.get_by_id(reaction)
            rxn.lower_bound, rxn.upper_bound = rbounds

    # Objective function and optimal percentage (this may overwrite bounds for objective function from previous dictionary)
    if objective_function is None:
        logging.warning("No objective function was provided. The first one as returned by the get_objective_functions() will be used.")
        objective_function = get_objective_functions(cobra_model)[0]

    modified_cobra_model.objective = {
        modified_cobra_model.reactions.get_by_id(objective_function): 1.0
    }
    modified_cobra_model.reactions.get_by_id(objective_function).upper_bound = 1000

    fba_solution = modified_cobra_model.optimize()
    objective    = fba_solution.objective_value

    if not 0 <= optimal_percentage <= 100:
        logging.error("Optimal percentage cannot be lower than 0 or greater than 100. Please provide a valid value.")
    fraction_of_optimum = float(optimal_percentage / 100)

    # In case numeric errors lead to a higher value
    minimum_objective_value = round((objective * fraction_of_optimum), 6)
    if minimum_objective_value > objective:
        minimum_objective_value = objective - 1e-3

    # Set minimum bound of objective value equal to the fraction of optimum provided
    modified_cobra_model.reactions.get_by_id(objective_function).lower_bound = minimum_objective_value

    # Set direction of the objective
    modified_cobra_model.objective_direction = objective_direction

    # Build a dingo model based on the modified one
    modified_dingo_model = MetabolicNetwork.from_cobra_model(modified_cobra_model)

    return modified_cobra_model, modified_dingo_model


def sample_optgp(cobra_model, n_samples = 3000, thinning = 100, reaction_in_rows = True):
    """
    Function that performs sampling with the OptGP sampler of cobrapy. Users can define number of samples,
    thinning parameter. Resulting samples dataset is returned in a numpy array with 
    reactions as rows or columns specified by user
    """

    sampler_optgp = OptGPSampler(cobra_model, thinning=thinning)
    samples_optgp = sampler_optgp.sample(n=n_samples)
    samples_optgp = samples_optgp.to_numpy()

    if reaction_in_rows == True:
        samples_optgp = samples_optgp.T

    return samples_optgp


def sample_dingo(dingo_model, reaction_in_rows = True, ess=1000, psrf = False, solver="gurobi", final_n_samples = None):
    """
    Function that performs sampling with dingo. Users can define the ess and psrf parameters.
    Resulting samples dataset is returned in a numpy array with reactions as rows or columns specified by user.
    User can call the "reduce_samples_dimensions" function to reduce the samples to a specified number
    """

    set_default_solver(solver)

    def reduce_samples_dimensions(samples, default_n_samples, final_n_samples):
        """
        Function that reduces the number of samples dingo returns to a user-specified value
        (current implementation may include the same sample more than 1 time in the final samples dataset)
        """

        if isinstance(final_n_samples, int):

            if final_n_samples < default_n_samples:

                # Reduce the dimensions of the samples to a specific value
                samples_initial_dimensions = samples.shape[1]

                step = samples_initial_dimensions / final_n_samples

                subsample_indices = np.arange(0, samples_initial_dimensions, step).astype(int)

                # force exactly the same number of samples. This may take duplicate samples
                subsample_indices = subsample_indices[:final_n_samples]

                final_samples = samples[:, subsample_indices]

            else:
                raise ValueError("samples_dimensions must be an integer and not exceed the current dimensions.")

        return final_samples

    sampler_dingo_default = PolytopeSampler(dingo_model)
    samples_dingo = sampler_dingo_default.generate_steady_states(ess=ess, psrf=psrf)
    default_n_samples = samples_dingo.shape[1]

    if final_n_samples != None:
        samples_dingo = reduce_samples_dimensions(samples_dingo, default_n_samples, final_n_samples)

    if reaction_in_rows == True:
        pass
    else:
        samples_dingo = samples_dingo.T

    return samples_dingo


def sample_gapsplit(cobra_model, n_samples = 3000, reaction_in_rows = True, add_loopless_cobrapy = False, fraction_of_optimum=0):
    """
    Function that performs sampling with the gapsplit algorithm. Users can define number of samples, if they want
    to modify their model with the "add_loopless" function and a value for the "fraction_of_optimum" parameter used in the gapsplit algorithm.
    Resulting samples dataset is returned in a numpy array with reactions as rows or columns specified by user
    """

    ec_cobra_model = cobra_model.copy()

    if add_loopless_cobrapy == True:
        add_loopless(ec_cobra_model)

    gapsplit_samples = gapsplit(
        ec_cobra_model,
        n=n_samples,
        gurobi_direct=False,
        fraction_of_optimum=fraction_of_optimum
    )

    gapsplit_samples = gapsplit_samples.to_numpy()

    if reaction_in_rows == True:
        gapsplit_samples = gapsplit_samples.T

    return gapsplit_samples


def export_to_pickle(samples, filename):
    """
    Function that exports samples to a pickle file
    """

    with open(filename, "wb") as samples_file:
        pickle.dump(samples, samples_file)


def load_from_pickle(filename):
    """
    Function that loads samples from a pickle file
    """

    with open(filename, "rb") as samples_file:
        samples_file = pickle.load(samples_file)

    return samples_file


def plot_grid_95_reactions(samples, cobra_model, tolerance=1e-3, nrows=20, ncols=5, title="Sampling Distributions"):
    """
    Function that plots a grid of the sampling distributions of all model reactions.
    User can define distributions plotted per rows and columns
    """

    cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]
    samples = samples.copy()

    if nrows * ncols < len(cobra_reactions_str):
        raise Exception("Change dimensions of the plot (nrows, ncols) to store all distributions")

    # if provided sampling dataset has reactions as rows ==> transpose
    if samples.shape[0] == len(cobra_reactions_str):
        samples = samples.T

    df = pd.DataFrame(samples)
    df.columns = cobra_reactions_str

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4 * ncols, 3 * nrows))
    axes = np.array(axes).flatten()

    for i, col in enumerate(df.columns):
        ax = axes[i]
        values = df[col].values
        range_width = values.max() - values.min()

        # Set transparency for very narrow distributions
        alpha = 0.4 if range_width < tolerance else 1.0

        ax.hist(values, bins=20, color='skyblue', edgecolor='black', alpha=alpha)
        ax.set_title(col)
        ax.grid(True)

    # Hide unused axes
    for j in range(len(df.columns), len(axes)):
        fig.delaxes(axes[j])

    fig.suptitle(title, fontsize=20, y=1.02)  # Title above all plots
    plt.tight_layout()
    plt.show()


def sampling_statistics(samples, reactions_ids_list=None, reaction_id=""):
    """
    Function that prints statistics for the sampling distribution of a specified model reaction
    """

    if reactions_ids_list == None:
        raise Exception("List with Reactions IDs not provided")

    if reaction_id in reactions_ids_list:
        reaction_index = reactions_ids_list.index(reaction_id)

        mean = np.mean(samples[reaction_index])
        min = np.min(samples[reaction_index])
        max = np.max(samples[reaction_index])
        std = np.std(samples[reaction_index])
        skewness = stats.skew(samples[reaction_index])
        kurtosis = stats.kurtosis(samples[reaction_index])

        return mean, min, max, std, skewness, kurtosis

    else:
        print("Reaction ID provided not in the list of Reactions IDs")
