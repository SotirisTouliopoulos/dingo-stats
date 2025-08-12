
# KEGG pathways information

## Map model reactions IDs to KEGG terms


The `map_model_to_kegg_reactions_dictionary` function will create a dictionary that will assign KEGG terms (values) to BiGG/SEED ids (keys) only from given model's information (without searching on online databases)

- `cobra_model` is a cobra model object

```python
initial_bigg_to_kegg_dictionary = map_model_to_kegg_reactions_dictionary(
    cobra_model = ec_cobra_model)

print(initial_bigg_to_kegg_dictionary.get("PFL"))
```

```python
'R00212'
```


The `read_json_file` function reads the JSON file where the KEGG information is stored and saves it in a JSON format and a pandas dataframe objects. This JSON file is available [here](https://raw.githubusercontent.com/MGXlab/DNNGIOR/refs/heads/main/docs/biochemistry/reactions.json) and has information that helps with mapping BiGG IDs, SEED IDs, KEGG terms and KEGG pathway names.

- `filepath` is the path where the `reactions.json` file is located

```python
reactions_json, reactions_pandas = read_json_file(
    filepath = "../ext_data/reactions/reactions.json")

print(reactions_pandas['aliases'][0])
print(reactions_pandas['linked_reaction'][0])
```

```python
['AraCyc: INORGPYROPHOSPHAT-RXN', 'BiGG: IPP1; PPA; PPA_1; PPAm', 'BrachyCyc: INORGPYROPHOSPHAT-RXN', 'KEGG: R00004', 'MetaCyc: INORGPYROPHOSPHAT-RXN', 'Name: Diphosphate phosphohydrolase; Inorganic diphosphatase; Inorganic pyrophosphatase; Pyrophosphate phosphohydrolase; diphosphate phosphohydrolase; inorganic diphosphatase; inorganic diphosphatase (one proton translocation); inorganicdiphosphatase; pyrophosphate phosphohydrolase']

rxn27946;rxn27947;rxn27948;rxn32487;rxn38157;rxn38158
```


The `dictionary_reaction_id_to_kegg_id` function given the pandas dataframe created from the `reactions.json` file, builds two dictionaries for fast lookup of KEGG reaction IDs from BiGG or SEED IDs. These dictionaries will be used as input in the `reaction_id_to_kegg_id` function.

- `reactions_pandas` is a pandas dataframe created form the `read_json_file` function and the `reactions.json` file

```python
bigg_to_kegg, seed_to_kegg = dictionary_reaction_id_to_kegg_id(
    reactions_pandas = reactions_pandas)

print(bigg_to_kegg.get("IPP1"))
print(seed_to_kegg.get("rxn27946"))
```

```python
R00004
R00004
```


The `reaction_id_to_kegg_id` function takes as arguments: a BiGG or a SEED id, the modeltype and the mapping dictionaries created above. It returns the corresponding KEGG id.

- `reaction_id` is the BiGG or SEED reaction ID (e.g., "SUCDi" or "rxn12345")
- `modeltype` is the model type (either "BiGG" or "SEED" and determines which dictionary to use)
- `bigg_to_kegg` is a dictionary mapping BiGG IDs to KEGG IDs
- `seed_to_kegg` is a dictionary mapping SEED IDs to KEGG IDs

```python
kegg_id_from_bigg = reaction_id_to_kegg_id(
    reaction_id = "IPP1", 
    modeltype = "BiGG", 
    bigg_to_kegg = bigg_to_kegg, 
    seed_to_kegg = seed_to_kegg)

kegg_id_from_seed = reaction_id_to_kegg_id(
    reaction_id = "rxn19264", 
    modeltype = "SEED", 
    bigg_to_kegg = bigg_to_kegg, 
    seed_to_kegg = seed_to_kegg)

print(kegg_id_from_bigg)
print(kegg_id_from_seed)
```

```python
R00004
R00009
```


The `fill_missing_kegg_ids_in_initial_dictionary` function is used to further map KEGG to BiGG/SEED ids, in cases where the initial model lacks some information. It fills in missing KEGG IDs (NAs) in the initial mapping dictionary

- `initial_bigg_to_kegg_dictionary` is a dictionary with reaction IDs as keys and KEGG IDs (or None) as values created from the `map_model_to_kegg_reactions_dictionary` function and includes only the default mapping information from the model
- `modeltype` is the model type (either "BiGG" or "SEED")
- `bigg_to_kegg` is a dictionary mapping BiGG IDs to KEGG IDs
- `seed_to_kegg` is a dictionary mapping SEED IDs to KEGG IDs

```python
final_bigg_to_kegg_dictionary = fill_missing_kegg_ids_in_initial_dictionary(
    initial_model_to_kegg_dictionary = initial_bigg_to_kegg_dictionary, 
    modeltype="BiGG",
    bigg_to_kegg = bigg_to_kegg,
    seed_to_kegg = seed_to_kegg)

print(initial_bigg_to_kegg_dictionary.get('PFK'))
print(final_bigg_to_kegg_dictionary.get('PFK'))
```

```python
None
R00756
```


The `get_kegg_pathways_from_reaction_ids` function fetches KEGG pathway information for a set of model reactions 
and creates a pandas dataframe with the following columns: with columns: 'model_reaction', 'kegg_reaction', 'pathway_ids', 'pathway_names'.

- `final_bigg_to_kegg_dictionary` is a dictionary with reaction IDs as keys and KEGG IDs (or None) as values created from the `fill_missing_kegg_ids_in_initial_dictionary` function and includes updated mapping information from the KEGG database.
- `max_workers` corresponds to the number of threads for parallel downloading

```python
df_kegg_pathways = get_kegg_pathways_from_reaction_ids(
    final_model_to_kegg_dictionary = final_bigg_to_kegg_dictionary,
    max_workers = 8)

print(df_kegg_pathways["model_reaction"].iloc[0])
print(df_kegg_pathways["kegg_reaction"].iloc[0])
print(df_kegg_pathways["pathway_ids"].iloc[0])
print(df_kegg_pathways["pathway_names"].iloc[0])
```

```python
PFL
R00212
[rn00620, rn00650, rn01100, rn01120]
[Pyruvate metabolism, Butanoate metabolism, ...]
```


## Subset reactions from pathways

The `subset_model_reactions_from_pathway_info` function given a dataFrame created wuth the `get_kegg_pathways_from_reaction_ids` function, returns all reaction IDs affiliated with a given KEGG pathway name or ID.

```python
PPP_from_name = subset_model_reactions_from_pathway_info(
    kegg_info_df = df_kegg_pathways, 
    pathway_query = "Pentose phosphate pathway")

Glycolysis_from_name = subset_model_reactions_from_pathway_info(
    kegg_info_df = df_kegg_pathways, 
    pathway_query = "Glycolysis / Gluconeogenesis")

Glycolysis_from_id = subset_model_reactions_from_pathway_info(
    kegg_info_df = df_kegg_pathways, 
    pathway_query = "rn00010")

print(PPP_from_name)
print(Glycolysis_from_name)
print(Glycolysis_from_id)
```

```python
['FBA', 'FBP', 'GND', 'PFK', 'PGL', 'RPE', 'RPI', 'TKT1']
['ALCD2x', 'ENO', 'FBA', 'FBP', 'GAPD', 'PFK', 'PGK', 'PGM', 'PPCK', 'PPS', 'PYK', 'TPI']
['ALCD2x', 'ENO', 'FBA', 'FBP', 'GAPD', 'PFK', 'PGK', 'PGM', 'PPCK', 'PPS', 'PYK', 'TPI']
```


The `sort_reactions_by_model_order` function flattens the lists provided in the `subsets` argument (corresponding to reactions from different pathways) in a single list and then orders the element of the new list based on the order of the reaction in the initial model. If any duplicates exist they are not removed by this function, so an additional step is required if user wants to exclude duplicate reaction IDs.

- `full_list` is the reference list that defines the desired order. Usually corresponds to the model reactions
- `*subsets` is/are one or more subset lists to be merged and ordered. Usually corresponds to reactions from pathways of interest.

```python
reactions_in_pathways_ordered_duplicates = sort_reactions_by_model_order(
    full_list = ec_cobra_reaction_ids, 
    Glycolysis,
    PPP)

# Additional step to remove duplicates
reactions_in_pathways_ordered = []
[reactions_in_pathways_ordered.append(val) for val in reactions_in_pathways_ordered_duplicates if val not in reactions_in_pathways_ordered]
```


The `dictionary_reaction_id_to_pathway` function takes one or multiple lists containing reaction IDs (corresponding to KEGG pathways and creates a dictionary that maps the IDs to pathway names. If a reaction appears in more than 1 pathway, it is classified with the term `Multiple-Pathways`. This is useful for plotting to work with subsets of reactions and to replace names from the `df_kegg_pathways` dataframe like `Glycolysis / Gluconeogenesis` to `Glycolysis` and `Pentose phosphate pathway` to `PPP`.

- `**named_lists` are named lists where each argument is a list of reaction IDs and the argument name represents the pathway name.

```python
bigg_to_pathway_dict = dictionary_reaction_id_to_pathway(
    Glycolysis = Glycolysis, 
    PPP = PPP)

print(bigg_to_pathway_dict.get("GND"))
print(bigg_to_pathway_dict.get("ENO"))
print(bigg_to_pathway_dict.get("FBA"))
```

```python
"Pentose phosphate pathway"
"Glycolysis / Gluconeogenesis"
"Multiple-Pathways"
```


The `reaction_in_pathway_binary_matrix` function is used to create a new pandas dataframe with reactions as rows and different pathways as columns. The corresponding cell of the dataframe will show if a reaction belongs to a certain pathway (1) or not (0). If a reaction belongs to more than one pathways, then the column: `Multiple-Pathways` is created and the reaction matching this will only get True (1) there and not in the individual pathway columns (e.g. 1 in `Multiple-Pathways`, 0 in `Glycolysis` and 0 in `PPP`).

- `reaction_id_to_pathway_dict` is dictionary mapping reaction IDs to pathway names created with the `dictionary_reaction_id_to_pathway` function

```python
binary_df = reaction_in_pathway_binary_matrix(
    reaction_id_to_pathway_dict = bigg_to_pathway_dict)
```


The `plot_reaction_in_pathway_heatmap` function is used to plot a heatmap of the `binary_df` created from the `reaction_in_pathway_binary_matrix` function to better illustrate the connection between reactions and pathways.

- `binary_df` is a pandas dataFrame with binary values (0 or 1)
- `font_size` is the font size for axis labels and ticks
- `fig_width` is the width of the figure in pixels
- `fig_height` is the height of the figure in pixels
- `title` is the title of the plot

```python
plot_reaction_in_pathway_heatmap(
    binary_df = binary_df, 
    font_size = 8, 
    fig_width = 600, 
    fig_height = 600, 
    title = "" )
```

![heatmap_pathways_binary](/img/heatmap_pathways_binary.png)


The `subset_sampling_array_from_reaction_ids` function subsets a sampling 2D array (with reactions as rows and samples as columns) to include only reactions of interest.

- `samples` is a sampling 2D array with reactions as rows and samples as columns
- `model_reactions` is a list containing the model's reactions
- `subset_reactions` is a list containing reactions of interest to subset the sampling array

```python
subset_pathways_optgp_condition_100 = subset_sampling_array_from_reaction_ids(
    samples = samples_optgp_condition_100, 
    model_reactions = ec_cobra_reaction_ids, 
    subset_reactions = reactions_in_pathways_ordered)

subset_pathways_optgp_condition_0 = subset_sampling_array_from_reaction_ids(
    samples = samples_optgp_condition_0, 
    model_reactions = ec_cobra_reaction_ids,
    subset_reactions = reactions_in_pathways_ordered)
```