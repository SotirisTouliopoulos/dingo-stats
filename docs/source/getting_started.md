
# Getting started

## Load, Modify and Sample models

The function `load_model` is used to load a metabolic model from a sbml format file and store it as a `cobra` model object
returned as first argument. Additionally the model's reactions objects and their corresponding string ids are returned
as second and third arguments.

```python
ec_cobra_model, ec_cobra_reactions, ec_cobra_reaction_ids,  = load_model("../ext_data/models/e_coli_core.xml")
```

The function `get_objective_functions` is used to get the assigned objective function(s) form the given cobra model

```python
objective_functions = get_objective_functions(ec_cobra_model)
```

The function `get_reaction_bounds` is used to get and save in a python dictionary the reaction ids as keys and the corresponding reaction bounds (lower/upper) as values.

```python
default_reaction_bounds = get_reaction_bounds(ec_cobra_model)
```

## Inspect Sampling distributions



