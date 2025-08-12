
# Getting started

## Load and Modify models


The function `load_model` is used to load a metabolic model (`cobra.Model` object) from a sbml format file and return it to the user (1st arg.) alongside the model's reactions (2nd arg.) and their corresponding string ids (3rd arg.).

- `filepath` is the path to the cobra model in sbml format

```python
ec_cobra_model, ec_cobra_reactions, ec_cobra_reaction_ids = load_model(
    filepath = "../ext_data/models/e_coli_core.xml")
```


The function `get_objective_functions` is used to get the assigned objective function(s) from the given cobra model

- `cobra_model` is the corresponding model in cobra format

```python
objective_functions = get_objective_functions(
    cobra_model = ec_cobra_model)

print(objective_functions)
```

```python
['BIOMASS_Ecoli_core_w_GAM']
```


The function `get_reaction_bounds` is used to create a python dictionary with reaction IDs as keys and the corresponding reaction bounds (lower/upper) as values.

- `cobra_model` is the corresponding model in cobra format

```python
default_reaction_bounds = get_reaction_bounds(
    cobra_model = ec_cobra_model)

print(default_reaction_bounds.get("BIOMASS_Ecoli_core_w_GAM"))
```

```python
(0.0, 1000.0)
```


The `modify_model` function modifies a given cobra model. Users can change objective function, optimal percentage (lower bound of the objective will be fixed based on the optimal percentage), define custom reaction bounds and change objective direction.

- `cobra_model` is a cobra model object
- `objective_function` is a string with the requested objective function ID to be assigned to the model
- `optimal_percentage` is the percentage of the FBA solution to be set as lower bound to the objectie function.
- `reaction_bounds` is a dictionary of custom reaction bounds to be assigned to the model's reactions
- `objective_direction` defines the direction to optimize the model: `max` and `min` are available options.

This function enables the user to create metabolic models simulating different conditions. Here, 2 new models are created based on updating the initial model:

- One that maximizes for biomass production (asking at least almost 100% of the biomass maximum FBA value)

```python
ec_cobra_model_condition_100 = modify_model(
    cobra_model         = ec_cobra_model,
    objective_function  = "BIOMASS_Ecoli_core_w_GAM",
    optimal_percentage  = 100,
    objective_direction = "max"
)
```

- One that maximize for biomass production (asking at least 0% of the biomass maximum FBA value)

```python
ec_cobra_model_condition_0 = modify_model(
    cobra_model         = ec_cobra_model,
    objective_function  = "BIOMASS_Ecoli_core_w_GAM",
    optimal_percentage  = 0,
    objective_direction = "max"
)
```


## Perform Sampling

For sampling, dingo-stats has various functions that integrate known samplers and make their outtput on a desired format for further analysis. However the user can use whichever sampler he wants, as long as final samples are concerted in a numpy 2D array.

The `sample_optgp` function performs sampling using the [OptGp](https://doi.org/10.1371/journal.pone.0086587) sampler. Users can change the `n_samples`, `thinning`, parameters:

- `cobra_model` is the corresponding model in cobra format
- `n_samples` -- defines the number of samples returned to the user.
- `thinning` -- defines the `thinning` parameter. Default to 100 means samples are returned every 100th step.
- `reaction_in_rows`  -- boolean variable that if `True` reactions are rows in the final sampling dataframe.

```python
samples_optgp_condition_100 = sample_optgp(
    cobra_model = ec_cobra_model_condition_100, 
    n_samples = 3000, 
    thinning=100, 
    reaction_in_rows = True)

samples_optgp_condition_0 = sample_optgp(
    cobra_model = ec_cobra_model_condition_0, 
    n_samples = 3000, 
    thinning=100, 
    reaction_in_rows = True)
```


The `sample_dingo` performs sampling with [dingo](https://doi.org/10.1093/bioadv/vbae037). Users can define the `ess` and `psrf` parameters.

- `dingo_model` is the corresponding model in dingo format
- `ess` stands for the effective sample size (ESS) (default value is 1000).
- `psrf` is a flag to request an upper bound equal to 1.1 for the value of the potential scale reduction factor of each marginal flux (default option is False).
- `solver` -- defines the linear programming solver
- `final_n_samples` -- integer that defines the final number of samples returned to the user (post sampling step)
- `reaction_in_rows` -- boolean variable that if `True` reactions are rows in the sampling dataframe.

First convert the modified cobra models to `dingo` models:

```python
from dingo import MetabolicNetwork

ec_dingo_model_condition_100 = MetabolicNetwork.from_cobra_model(
    ec_cobra_model_condition_100)

ec_dingo_model_condition_0 = MetabolicNetwork.from_cobra_model(
    ec_cobra_model_condition_0)
```

And then perform sampling on the dingo models:

```python
samples_dingo_condition_100 = sample_dingo(
    dingo_model = ec_dingo_model_condition_100, 
    reaction_in_rows = True, 
    ess=1000, 
    psrf = False, 
    solver="gurobi", 
    final_n_samples = None)

samples_dingo_condition_0 =   sample_dingo(
    dingo_model = ec_dingo_model_condition_0, 
    reaction_in_rows = True, 
    ess=1000, 
    psrf = False, 
    solver="gurobi", 
    final_n_samples = None)
```


## Inspect Sampling distributions

The `plot_grid_reactions_flux` function plots a grid of selected flux distributions and enables inspections of sampling results to get early insights. Distributions here are chosen based on a list of IDs derived from KEGG pathway information. The corresponding KEGG pathways tutorial is presented in the next section. From the grid we can detect: Left/Right Shifted, Normal, Fixed (transparent based on the `tolerance` parameter) or other uncommon distributions.

- `samples` is a numpy 2D array of the samples
- `model_reactions` is a list containing strings of reaction IDs
- `tolerance` is a tolerance level to make distribution plot transparent based on flux range
- `nrows` is the number of rows for the plot
- `ncols` is the number of columns for the plot
- `title` is the title of the plot

In the following chunk the grid from Glycolysis and Pentose-Phosphate Pathway (PPP) pathways is shown:

```python
plot_grid_reactions_flux(
    samples=samples_dingo_condition_100,
    model_reactions=ec_cobra_reaction_ids, 
    nrows=5, 
    ncols=5,
    tolerance=1e-3,
    title="Sampling Distributions")
```

![grid_flux_distributions](/img/grid_flux_distributions.png)


The `sampling_statistics` function calculate statistics on flux distributions for a reaction of interest. This information can be used to identify significantly altered flux distributions between different conditions

- `samples` is a numpy 2D array of the samples
- `model_reactions` is a list containing strings of reaction IDs
- `reaction_id` is a reaction ID of the selected reaction to calculate statistics on

```python
mean, min, max, std, skewness, kurtosis = sampling_statistics(
    samples = samples_optgp_condition_100, 
    model_reactions=ec_cobra_reaction_ids,
    reaction_id="FRD7")

print(mean, min, max, std, skewness, kurtosis)
```

```python
479.1369619828837 0.14995320273677437 993.9157983127856 289.3315287517938 0.09951034266820558 -1.2064534542990524
```

```python
mean, min, max, std, skewness, kurtosis = sampling_statistics(
    samples = samples_optgp_condition_0,
    model_reactions=ec_cobra_reaction_ids,
    reaction_id="FRD7")

print(mean, min, max, std, skewness, kurtosis)
```

```python
474.7669865820229 0.10135434309104707 997.5267773446074 278.80950212031655 0.07570051318545512 -1.1733598969449195
```