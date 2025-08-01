
from statsmodels.stats.multitest import multipletests
from scipy import stats
import numpy as np
from typing import Dict, Tuple, List
from cobra import Model, Reaction
from numpy.typing import NDArray
import pandas as pd
from scipy.stats import hypergeom
import plotly.express as px


def significantly_altered_reactions(
    conditions: List[NDArray[np.float64]] = [], 
    selected_comparisons: Tuple = [(0, 1)], 
    cobra_model: Model = None,
    p_value_cutoff: float = 0.05,
    fold_change_cutoff: float = 1.2,
    std_cutoff: float = 1e-2
) -> Tuple[List[str], List[str]]:

    """
    Function that takes as input at least 2 flux sampling conditions to compare and identifies significantly altered reactions
    It performs a Kolmogorov-Smirnov (KS) non-parametric test and corrects p-value for multiple comparisons.
    It additionally calculates a fold change that together with the p-value classifies reactions as significantly altered or not. 
    
    Keyword arguments:
    conditions (List) -- List of different flux sampling arrays
    selected_comparisons (List) -- List showing which conditions to compare (useful when comparing more than 2 sampling arrays)
    cobra_model (Model) -- cobra model object
    p_value_cutoff (float) -- cutoff for p value from KS test to consider 2 distributions significantly different
    fold_change_cutoff (float) -- cutoff for fold-change to consider 2 distributions significantly different
    std_cutoff (float) -- cutoff to ensure distributions that are compared are not fixed to a certain value

    Returns:
    Tuple[List, List]
        significant_diff_reactions (List) -- List containing reactions significantly altered between the conditions
        not_significant_diff_reactions (List) -- List containing reactions not significantly altered between the conditions
        pval_dict (Dict) -- Dictionary mapping reaction ID to corrected p-value.
        fold_change_dict (Dict) -- Dictionary mapping reaction ID to fold-change value.
    """

    cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]
        
    # if provided sampling dataset has reactions as cols ==> transpose
    for condition in conditions:
        if condition.shape[0] == condition.shape[1]:
            raise ValueError("Samples array provided has equal rows and columns dimensions. Please change the number of samples")
        
        if condition.shape[1] == len(cobra_reactions_str):
            condition = condition.T

    # Store p-values for each row for KS test
    p_values = {row: [] for row in range(len(cobra_reactions_str)) }
    # Store ks-values for each row for KS test
    ks_values = {row: [] for row in range(len(cobra_reactions_str)) }
    # Store fold-change-values for each row for KS test
    fold_change_values = {row: [] for row in range(len(cobra_reactions_str)) }

    for row in range(len(cobra_reactions_str)):
        for i, j in selected_comparisons:
            
            fold_change = np.absolute( (np.mean(conditions[i][row]) - np.mean(conditions[j][row]) ) / (np.mean(conditions[j][row]) + 1e-8) )
            #fold_change = (np.mean(conditions[i][row]) - np.mean(conditions[j][row]) ) / (np.mean(conditions[j][row]) + 1e-8)
            fold_change_values[row].append(fold_change)

            if (np.std(conditions[i][row]) > std_cutoff) or (np.std(conditions[j][row]) > std_cutoff):

                ks, p = stats.ks_2samp(conditions[i][row], conditions[j][row], alternative='two-sided')
                p_values[row].append(p)
                ks_values[row].append(ks)

            else:
                p_values[row].append(1)
                ks_values[row].append(1)


    p_values_copy = list(p_values.values())
    flat_p_values = np.array(p_values_copy).flatten()
    # Apply FDR correction
    _, corrected_p_values, _, _ = multipletests(flat_p_values, method='fdr_bh')
    # Reshape corrected p-values back to the original matrix shape
    corrected_p_values = corrected_p_values.reshape(np.array(p_values_copy).shape)
    
    fold_change_values = np.array(list(fold_change_values.values()))

    significant_diff_indices = np.where(
    np.logical_and(
        corrected_p_values < p_value_cutoff,
        np.abs(fold_change_values) > fold_change_cutoff
    ) )[0]

    significant_diff_reactions = [cobra_reactions_str[i] for i in significant_diff_indices]

    not_significant_diff_reactions = list(set(cobra_reactions_str) - set(significant_diff_reactions))
    
    pval_dict = {cobra_reactions_str[i]: corrected_p_values[i][0] for i in range(len(cobra_reactions_str))}
    fold_change_dict = {cobra_reactions_str[i]: fold_change_values[i][0] for i in range(len(cobra_reactions_str))}

    return significant_diff_reactions, not_significant_diff_reactions, pval_dict, fold_change_dict


def plot_volcano(
    pval_dict: dict,
    fold_change_dict: dict,
    p_value_cutoff: float = 0.05,
    fold_change_cutoff: float = 1.2,
    annotate: list = [],
    width: int = 800,
    height: int = 600,
    title: str = "",
    show_cutoff_lines: bool = True
):
    """
    Creates a volcano plot to show results from differential flux analysis.

    Keyword arguments:
        pval_dict (dict) -- Dictionary mapping reaction ID to corrected p-value.
        fold_change_dict (dict) -- Dictionary mapping reaction ID to fold-change value.
        p_value_cutoff (float) -- Significance threshold for p-value.
        fold_change_cutoff (float) -- Threshold for fold-change.
        annotate (list) -- List of reaction IDs to annotate on the plot.
        width (int) -- Width of the figure in pixels.
        height (int) -- Height of the figure in pixels.
        title (str) -- Title of the plot.
        show_cutoff_lines (bool) -- Whether to draw p-value and fold-change cutoff lines.
    """
    
    df = pd.DataFrame({
        'reaction': list(pval_dict.keys()),
        'p_value': [pval_dict[r] for r in pval_dict],
        'fold_change': [fold_change_dict[r] for r in pval_dict],
    })

    df['-log10_pval'] = -np.log10(df['p_value'] + 1e-10)
    
    def significance(row):
        if row['p_value'] < p_value_cutoff and abs(row['fold_change']) > fold_change_cutoff:
            return 'Significant'
        return 'Not Significant'
    
    df['significance'] = df.apply(significance, axis=1)

    fig = px.scatter(
        df,
        x='fold_change',
        y='-log10_pval',
        color='significance',
        hover_data=['reaction'],
        title=title,
        color_discrete_map={'Significant': 'red', 'Not Significant': 'gray'},
        width=width,
        height=height
    )

    # Annotate reactions if specified
    if annotate:
        for reaction in annotate:
            if reaction in df['reaction'].values:
                point = df[df['reaction'] == reaction].iloc[0]
                fig.add_annotation(
                    x=point['fold_change'],
                    y=point['-log10_pval'],
                    text=reaction,
                    showarrow=True,
                    arrowhead=1,
                    ax=20,
                    ay=-20
                )

    # Add cutoff lines if requested
    if show_cutoff_lines:
        fig.add_shape(
            type='line',
            x0=np.log2(fold_change_cutoff), x1=np.log2(fold_change_cutoff),
            y0=0, y1=df['-log10_pval'].max(),
            line=dict(dash='dash', color='blue'),
        )
        fig.add_shape(
            type='line',
            x0=-np.log2(fold_change_cutoff), x1=-np.log2(fold_change_cutoff),
            y0=0, y1=df['-log10_pval'].max(),
            line=dict(dash='dash', color='blue'),
        )
        fig.add_shape(
            type='line',
            x0=df['fold_change'].min(), x1=df['fold_change'].max(),
            y0=-np.log10(p_value_cutoff), y1=-np.log10(p_value_cutoff),
            line=dict(dash='dash', color='green'),
        )

    fig.update_layout(
        xaxis_title='Fold Change',
        yaxis_title='-log10(p-value)',
        template='plotly_white'
    )
    
    fig.show()


def dictionary_reaction_to_all_pathways(
    df_kegg_pathways: pd.DataFrame,
    reaction_col: str = 'model_reaction', 
    pathway_col: str = 'pathway_names'
) -> Dict:
    """
    Function that builds a dictionary mapping each reaction to a list of all the pathway names it belongs.

    Keyword arguments:
    df_kegg_pathways (pd.DataFrame): pandas DataFrame with reactions and pathway names.
    reaction_col: column name for reaction identifiers.
    pathway_col: column name for list of pathways.

    Returns:
    reaction_to_pathways (dict): Dictionary with reaction IDs as keys and pathway names as values: {reaction_id: [pathway1, pathway2, ...]}
    """

    reaction_to_pathways = {}
    
    for _, row in df_kegg_pathways.iterrows():
        reaction = row[reaction_col]
        pathways = row[pathway_col]
        
        # Ensure pathways is a list and not empty
        if isinstance(pathways, list) and pathways:
            reaction_to_pathways[reaction] = pathways

    return reaction_to_pathways


def hypergeometric_test_pathway_enrichment(
    significant_reactions: List,
    model_reactions: List, 
    reaction_to_pathways: Dict,
) -> pd.DataFrame:
    """
    Function that perform a hypergeometric test to find significantly affected pathways between our sampling conditions.

    Keyword arguments:
    significant_reactions (List): List of reaction IDs with altered flux.
    model_reactions (List) -- List of all reaction IDs considered in the analysis.
    reaction_to_pathways (Dict): Dictinary mapping reaction ID -> list of pathway names.

    Returns:
    hypergeometric_enrichment_df (pd.DataFrame) -- DataFrame with pathway, p-value, counts and adjusted FDR.
    """
    
    def add_enrichment_metrics(hypergeometric_pathway_enrichment_df: pd.DataFrame) -> pd.DataFrame:
        """
        Helper function to add fold enrichment and -log10(p-value) columns to the enrichment results DataFrame.
        
        Keyword arguments:
        hypergeometric_pathway_enrichment_df (pd.DataFrame) -- DataFrame from hypergeometric_pathway_enrichment
        
        Returns:
        hypergeometric_pathway_enrichment_df (pd.DataFrame) -- Updated DataFrame with ofld enrichment and -log10(p-value) information
        """
        
        #enrichment_metrics_df = hypergeometric_pathway_enrichment_df.copy()
        hypergeometric_pathway_enrichment_df['fold_enrichment'] = (hypergeometric_pathway_enrichment_df['significant_reactions_in_pathway'] / hypergeometric_pathway_enrichment_df['reactions_in_pathway']) / (hypergeometric_pathway_enrichment_df['all_significant_reactions'] / hypergeometric_pathway_enrichment_df['all_reactions'])
                                                
        hypergeometric_pathway_enrichment_df['log10_pval'] = -np.log10(hypergeometric_pathway_enrichment_df['pval'].clip(lower=1e-300))
        return hypergeometric_pathway_enrichment_df
    

    # Convert inputs to sets for efficiency
    significant_reactions = set(significant_reactions)
    model_reactions = set(model_reactions)

    # Build pathway -> reactions mapping
    pathway_to_reactions = {}
    for rxn, pathways in reaction_to_pathways.items():
        for pw in pathways:
            pathway_to_reactions.setdefault(pw, set()).add(rxn)

    results = []
    # Total number of reactions
    M = len(model_reactions)
    # Number of significant reactions
    n = len(significant_reactions)

    for pathway, pathway_reactions in pathway_to_reactions.items():
        # Reactions in this pathway
        K = len(model_reactions & pathway_reactions)
        # Significant reactions in this pathway
        k = len(significant_reactions & pathway_reactions)

        # Skip if the pathway has no reactions in the background
        if K == 0:
            continue

        # Hypergeometric test
        pval = hypergeom.sf(k - 1, M, n, K)

        results.append({
            'pathway': pathway,
            'significant_reactions_in_pathway': k,
            'reactions_in_pathway': K,
            'all_significant_reactions': n,
            'all_reactions': M,
            'pval': pval
        })

    # Multiple testing correction (Benjamini-Hochberg FDR)
    hypergeometric_enrichment_df = pd.DataFrame(results)
    if not hypergeometric_enrichment_df.empty:
        hypergeometric_enrichment_df['fdr'] = multipletests(hypergeometric_enrichment_df['pval'], method='fdr_bh')[1]
        hypergeometric_enrichment_df.sort_values('pval', inplace=True)
        

    hypergeometric_enrichment_df = add_enrichment_metrics(hypergeometric_enrichment_df)

    return hypergeometric_enrichment_df


def plot_pathway_enrichment(
    hypergeometric_enrichment_df: pd.DataFrame, 
    pval_threshold: float = 0.05, 
    use_fdr: bool = True,
    font_size: int = 14,
    width: int = 900,
    height: int = 600
):
    """
    Function that takes as input the `hypergeometric_pathway_enrichment_df` dataframe created with the 
    `hypergeometric_test_pathway_enrichment` function and plots the enrichment results in a bubble plot. 
    Bubble size stands for pathway size (number of reactions) and reaction colour stands for -log10(p-value)

    Keyword arguments:
    hypergeometric_enrichment_df (pd.DataFrame): Enrichment DataFrame from the `hypergeometric_test_pathway_enrichment` function.
    pval_threshold (float): significance p value threshold for filtering.
    use_fdr (bool): boolean for whether to use 'fdr' or 'pval' for filtering.
    font_size (int): font size for plot labels and title.
    width (int): width of the figure in pixels.
    height (int): height of the figure in pixels.
    """

    df = hypergeometric_enrichment_df.copy()

    # Filter based on significance
    pval_col = 'fdr' if use_fdr else 'pval'
    df = df[df[pval_col] < pval_threshold]

    # Sort for better visual order
    df = df.sort_values('fold_enrichment', ascending=False)

    fig = px.scatter(
        df,
        x='fold_enrichment',
        y='pathway',
        size='significant_reactions_in_pathway',
        color='log10_pval',
        color_continuous_scale='Reds',
        size_max=20,
        labels={
            'fold_enrichment': 'Fold Enrichment',
            'log10_pval': '-log10(p-value)',
            'significant_reactions_in_pathway': '# Reactions',
            'pathway': 'Pathway'
        },
        title='Pathway Enrichment Bubble Plot'
    )

    fig.update_layout(
        yaxis=dict(autorange="reversed"),
        font=dict(size=font_size),
        width=width,
        height=height
    )

    fig.show()