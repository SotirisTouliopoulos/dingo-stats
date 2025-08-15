
# Differential Flux Analysis

## Compare reaction flux distributions across conditions with a Kolmogorov–Smirnov test

The `significantly_altered_reactions` function takes as input at least 2 flux sampling conditions to compare and identifies significantly altered reactions. It performs a Kolmogorov-Smirnov (KS) non-parametric test and corrects p-value for multiple comparisons. It additionally calculates a fold change that together with the p-value classifies reactions as significantly altered or not. It returns 2 lists one with the significantly changed (`significant_diff_reactions`) and one with the not significantly changed (`not_significant_diff_reactions`) reactions. Also, 2 dicitonaries mapping reaction IDs to corrected p_values (`pval_dict`) and fold change values (`fold_change_dict`)

- `conditions` is a list of different flux sampling arrays
- `selected_comparisons` is a list showing which conditions to compare (useful when comparing more than 2 sampling arrays)
- `cobra_model` is a cobra model object
- `fold_change_cutoff` is a cutoff for fold-change to consider 2 distributions significantly different
- `std_cutoff` is a cutoff to ensure distributions that are compared are not fixed to a certain value

```python
conditions = [samples_optgp_condition_100, samples_optgp_condition_0]

selected_comparisons = [(0, 1)]

(significant_diff_reactions,
 not_significant_diff_reactions,
 pval_dict,
 fold_change_dict) = significantly_altered_reactions(
    conditions = conditions, 
    selected_comparisons = selected_comparisons, 
    cobra_model = ec_cobra_model,
    fold_change_cutoff = 0.6,
    std_cutoff = 1e-3)
```


The `significantly_altered_reactions` function creates a volcano plot to show results from differential flux analysis. Volcano plot has Fold Change on the x-axis and -log10(p_value) on the y-axis. Users can provide a reactions list in the `annotate` parameter, to show reaction IDs on the plot. Also, lines showing the significance cutoffs may optionally be added, when providing `p_value_cutoff`, `fold_change_cutoff` and `show_cutoff_lines`.

- `pval_dict` is a dictionary mapping reaction ID to corrected p-value.
- `fold_change_dict` is a dictionary mapping reaction ID to fold-change value.
- `p_value_cutoff` is a significance threshold for p-value (to appear in the plot).
- `fold_change_cutoff` is a threshold for fold-change (to appear in the plot).
- `annotate` is a list of reaction IDs to annotate on the plot.
- `width` is the width of the figure in pixels.
- `height` is the height of the figure in pixels.
- `title` is the title of the plot.
- `show_cutoff_lines` is a boolean variable for whether to draw p-value and fold-change cutoff lines.

```python
reactions_to_annotate = ["NH4t"]

plot_volcano(
    pval_dict = pval_dict,
    fold_change_dict = fold_change_dict,
    p_value_cutoff = 0.05,
    fold_change_cutoff = 0.6,
    annotate = reactions_to_annotate,
    width = 800,
    height = 600,
    title = "",
    show_cutoff_lines = True)
```

![volcano_plot](/img/volcano_plot.png)



## Hypergeometric test for pathway enrichment

The `dictionary_reaction_to_all_pathways` function is used to create a dictionary mapping reaction IDs to all pathways they belong.

- `df_kegg_pathways` is a pandas DataFrame with reactions and pathway names, created from the `get_kegg_pathways_from_reaction_ids` function
- `reaction_col` is the column name for reaction identifiers, as expected to appear in the `df_kegg_pathways` dataframe
- `pathway_col` is the column name for list of pathways, as expected to appear in the `df_kegg_pathways` dataframe

```python
reaction_to_pathways = dictionary_reaction_to_all_pathways(
    df_kegg_pathways = df_kegg_pathways,
    reaction_col = 'model_reaction', 
    pathway_col = 'pathway_names')

print(reaction_to_pathways.get("PGI"))
```

```python
['Starch and sucrose metabolism', 'Metabolic pathways', 'Biosynthesis of secondary metabolites']
```


The `hypergeometric_test_pathway_enrichment` function performs a hypergeometric test to find significantly affected pathways between our sampling conditions. It also calculated fold enrichment and p_values useful for filtering significant changes

- `significant_reactions` is a list of reaction IDs with altered flux.
- `model_reactions` is a list of all reaction IDs considered in the analysis.
- `reaction_to_pathways` is a dictinary mapping reaction ID -> list of pathway names.

```python
hypergeometric_pathway_enrichment_df = hypergeometric_test_pathway_enrichment(
    significant_reactions = significant_diff_reactions, 
    model_reactions = ec_cobra_reaction_ids, 
    reaction_to_pathways = reaction_to_pathways)
```


The `plot_pathway_enrichment` function takes as input the `hypergeometric_pathway_enrichment_df` dataframe created with the `hypergeometric_test_pathway_enrichment` function and plots the enrichment results in a bubble plot. Bubble size stands for pathway size (number of reactions) and reaction colour stands for -log10(p-value)

- `hypergeometric_enrichment_df` is a dataframe from the `hypergeometric_test_pathway_enrichment` function.
- `pval_threshold` is a significance p value threshold for filtering.
- `use_fdr` is a boolean for whether to use 'fdr' or 'pval' for filtering.
- `font_size` is the font size for plot labels and title.
- `width` is the width of the figure in pixels.
- `height` is the height of the figure in pixels.

```python
plot_pathway_enrichment(
    hypergeometric_enrichment_df = hypergeometric_pathway_enrichment_df, 
    pval_threshold = 2,
    use_fdr = True,
    font_size = 14,
    width = 900,
    height = 600)
```

![pathway_enrichment](/img/pathway_enrichment.png)
