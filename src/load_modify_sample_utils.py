
from dingo import MetabolicNetwork, PolytopeSampler
from dingo import set_default_solver
from cobra.io import read_sbml_model
from cobra.sampling.optgp import OptGPSampler
from cobra.util.solver import linear_reaction_coefficients
from gapsplit import gapsplit
from cobra.flux_analysis.loopless import add_loopless
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import pickle
import logging
from typing import Dict, Tuple, List
from cobra import Model, Reaction
from numpy.typing import NDArray


def load_model(filepath: str) -> Tuple[Model, List[Reaction], List[str]]:
    """
    Function to load and return a model in cobra format from given path.
    It returns 3 objects, the cobra model and 2 lists, one with the model reactions (cobra_reactions)
    and one with their IDs (cobra_reactions_ids)

    Keyword arguments:
    filepath -- path to cobra model in sbml format
    
    Returns:
    Tuple[Model, List, List]:
        cobra_model -- metabolic model in cobra format
        cobra_reactions -- list with the model reactions
        cobra_reactions_ids -- list with the string IDs of model reactions
    """

    cobra_model = read_sbml_model(filepath)
    cobra_reactions = cobra_model.reactions
    cobra_reactions_ids = [reaction.id for reaction in cobra_reactions]

    return cobra_model, cobra_reactions, cobra_reactions_ids


def get_objective_functions(cobra_model: Model) -> List[str]:
    """
    Function to find and return the objective functions of a cobra model
    
    Keyword arguments:
    cobra_model (Model) -- cobra model object
    
    Returns:
    objective_functions_ids (List) -- List containing one or multiple objective functions of given model
    """

    objectives_dict         = linear_reaction_coefficients(cobra_model)
    objective_functions     = list(objectives_dict.keys())
    objective_functions_ids = [rxn.id for rxn in objective_functions]

    return objective_functions_ids


def get_reaction_bounds(cobra_model: Model) -> Dict[str, Tuple[float, float]]:
    """
    Function to create and return a dictionary of a cobra's model reaction bounds.
    Reaction IDs are keys and their corresponding bounds are values.
    
    Keyword arguments:
    cobra_model (Model) -- cobra model object
    
    Returns:
    reaction_bounds_dict (Dict) -- Dictionary containing reaction IDs as keys and their bounds as values
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
    cobra_model: Model, 
    objective_function: str = None,
    optimal_percentage: float = 100,
    reaction_bounds: dict = {},
    objective_direction: str = "max"
) -> Model:

    """
    Function that modifies a given cobra model.

    Users can change objective function, optimal percentage, define custom reaction bounds and change objective 
    direction. Lower bound of the objective will be fixed based on the optimal percentage, in all cases. 
    
    Keyword arguments:
    cobra_model (Model) -- cobra model object
    objective_function (str) -- string with the requested objective function ID to be assigned to the model
    optimal_percentage (float) -- percentage of the FBA solution to be set as lower bound to the objectie function.
    reaction_bounds (dict) -- dictionary of custom reaction bounds to be assigned to the model's reactions
    objective_direction (str) -- the direction to optimize the model: `max` and `min` are available options.
    
    Returns:
    modified_cobra_model (Model) -- A modified cobra model object
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

    # in case numeric errors lead to a fraction_of_optimum value higher than 1
    if fraction_of_optimum > 1:
        fraction_of_optimum = 1

    minimum_objective_value = round( (objective * fraction_of_optimum), 6)
    if minimum_objective_value > objective:
        minimum_objective_value = objective - 1e-3

    modified_cobra_model.reactions.get_by_id(objective_function).lower_bound = minimum_objective_value
    modified_cobra_model.objective_direction = objective_direction

    return modified_cobra_model


def sample_optgp(cobra_model: Model,
                 n_samples: int = 3000,
                 thinning: int = 100,
                 reaction_in_rows: bool = True
) -> NDArray[np.float64]:
    """
    Function that performs sampling with the OptGP sampler of cobrapy. Users can define number of samples,
    thinning parameter. Resulting samples dataset is returned in a numpy array with 
    reactions as rows or columns specified by user
    
    Keyword arguments:
    cobra_model (Model) -- cobra model object
    n_samples (int) -- defines the number of samples returned to the user.
    thinning (int) -- defines the `thinning` parameter. Default to 100 means samples are returned every 100th step.
    reaction_in_rows (bool) -- boolean variable that if `True` reactions are rows in the sampling dataframe.
    
    Returns:
    samples_optgp (NDArray[np.float64]) -- Numpy 2D array of the sampling results
    """

    sampler_optgp = OptGPSampler(cobra_model, thinning=thinning)
    samples_optgp = sampler_optgp.sample(n=n_samples)
    samples_optgp = samples_optgp.to_numpy()

    if reaction_in_rows == True:
        samples_optgp = samples_optgp.T

    return samples_optgp


def sample_dingo(dingo_model, 
                 ess: int = 1000,
                 psrf: bool = False,
                 solver: str = "gurobi",
                 final_n_samples: int = None,
                 reaction_in_rows: bool = True
) -> NDArray[np.float64]:
    """
    Function that performs sampling with dingo. Users can define the ess and psrf parameters.
    Resulting samples dataset is returned in a numpy array with reactions as rows or columns specified by user.
    User can call the "reduce_samples_dimensions" function to reduce the samples to a specified number
    
    Keyword arguments:
    dingo_model -- dingo model object
    ess -- stands for the effective sample size (ESS) (default value is 1000).
    psrf -- is a flag to request an upper bound equal to 1.1 for the value of the potential scale reduction factor of each marginal flux (default option is False).
    solver -- defines the linear programming solver
    final_n_samples -- integer that defines the final number of samples returned to the user (post sampling step)
    reaction_in_rows -- boolean variable that if `True` reactions are rows in the sampling dataframe.
    
    Returns:
    samples_dingo (NDArray[np.float64]) -- Numpy 2D array of the sampling results
    """

    set_default_solver(solver)


    def reduce_samples_dimensions(samples: NDArray[np.float64], 
                                  default_n_samples: int, 
                                  final_n_samples: int
    ) -> NDArray[np.float64]:
        """
        Function that reduces the number of samples dingo returns to a user-specified value
        (current implementation may include the same sample more than 1 time in the final samples dataset)
        
        Keyword arguments:
        samples (NDArray[np.float64]) -- initial samples to be subset
        default_n_samples (int) -- number of initial samples
        final_n_samples (int) -- number of desired final samples
        
        Returns:
        final_samples (NDArray[np.float64]) -- Numpy 2D array of the subset sampling results
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


def sample_gapsplit(cobra_model: Model,
                    n_samples: int = 3000,
                    reaction_in_rows: bool = True,
                    add_loopless_cobrapy: bool = False,
                    fraction_of_optimum: float = 0
) -> NDArray[np.float64]:
    """
    Function that performs sampling with the gapsplit algorithm. Users can define number of samples, if they want
    to modify their model with the "add_loopless" function and a value for the "fraction_of_optimum" parameter used in the gapsplit algorithm.
    Resulting samples dataset is returned in a numpy array with reactions as rows or columns specified by user
    
    Keyword arguments:
    cobra_model (Model) -- cobra model object
    n_samples (int) -- defines the number of samples returned to the user.
    reaction_in_rows (bool) -- boolean variable that if `True` reactions are rows in the sampling dataframe
    add_loopless_cobrapy (bool) -- Boolean variable that if True modifies the model with the `add_loopless` function
    fraction_of_optimum (float) -- Float variable that defines the `fraction_of_optimum` parameter used in `gapsplit` algorithm
    
    Returns:
    gapsplit_samples (NDArray[np.float64]) -- Numpy 2D array of the sampling results
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


def export_to_pickle(samples: NDArray[np.float64], filename: str):
    """
    Function that exports samples to a pickle file
    
    Keyword arguments: 
    samples (NDArray[np.float64]) -- Numpy 2D array of the samples
    filename (str) -- String to name the created file
    """

    with open(filename, "wb") as samples_file:
        pickle.dump(samples, samples_file)


def load_from_pickle(filename: str) -> NDArray[np.float64]:
    """
    Function that loads samples from a pickle file
    
    Keyword arguments: 
    filename (str) -- String of the filename containing the samples
    """

    with open(filename, "rb") as samples_file:
        samples_file = pickle.load(samples_file)

    return samples_file


def plot_grid_reactions_flux(samples: NDArray[np.float64],
                             model_reactions: List[str],
                             tolerance: float = 1e-3,
                             nrows: int = 20,
                             ncols: int = 5,
                             title="Sampling Distributions"):
    """
    Function that plots a grid of the sampling distributions of all model reactions.
    User can define distributions plotted per rows and columns
    
    Keyword arguments:
    samples (NDArray[np.float64]) -- Numpy 2D array of the samples
    model_reactions (List[str]) -- A list containing strings of reaction IDs
    tolerance (float) -- tolerance level to make distribution plot transparent based on flux range
    nrows (int) -- Number of rows for the plot
    ncols (int) -- Number of columns for the plot
    title (str) -- Title of the plot
    """

    samples = samples.copy()

    if nrows * ncols < len(model_reactions):
        raise Exception("Change dimensions of the plot (nrows, ncols) to store all distributions")

    # if provided sampling dataset has reactions as rows ==> transpose
    if samples.shape[0] == len(model_reactions):
        samples = samples.T

    df = pd.DataFrame(samples)
    df.columns = model_reactions

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


def sampling_statistics(samples: NDArray[np.float64],
                        model_reactions: List[str] = None, 
                        reaction_id: str = ""
) -> Tuple[float, float, float, float, float, float]:
    """
    Function that prints statistics for the sampling distribution of a specified model reaction
    
    Keyword arguments:
    samples (NDArray[np.float64]) -- Numpy 2D array of the samples
    model_reactions (List[str]) -- A list containing strings of reaction IDs
    reaction_id (str) -- Reaction ID to calculate sampling statistics on
    
    Returns:
    Tuple[float, float, float, float, float, float]
        mean (float) -- average value from the reaction's flux distribution
        min (float) -- minimum value from the reaction's flux distribution
        max (float) -- maximum value from the reaction's flux distribution
        std (float) -- standard deviation value from the reaction's flux distribution
        skewness (float) -- skewness value from the reaction's flux distribution
        kurtosis (float) -- kurtosis value from the reaction's flux distribution
    """

    if model_reactions == None:
        raise Exception("List with Reactions IDs not provided")

    if reaction_id in model_reactions:
        reaction_index = model_reactions.index(reaction_id)

        mean = np.mean(samples[reaction_index])
        min = np.min(samples[reaction_index])
        max = np.max(samples[reaction_index])
        std = np.std(samples[reaction_index])
        skewness = stats.skew(samples[reaction_index])
        kurtosis = stats.kurtosis(samples[reaction_index])

        return mean, min, max, std, skewness, kurtosis

    else:
        print("Reaction ID provided not in the list of Reactions IDs")