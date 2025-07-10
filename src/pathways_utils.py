
import json
import pandas as pd
from bioservices import KEGG
import plotly.graph_objects as go


def read_json_file(filepath):
    with open(filepath, 'r') as file:
        reactions_json = json.load(file)
        reactions_pandas = pd.DataFrame(reactions_json)

        return reactions_json, reactions_pandas
    

def map_model_to_kegg_reactions(model):
    bigg_to_kegg_dictionary = {
        reac.id: reac.annotation.get("kegg.reaction") for reac in model.reactions}

    return bigg_to_kegg_dictionary


def bigg_to_kegg_id(bigg_id, df):
    """
    Given a BiGG ID, return the corresponding KEGG ID (from abbreviation column) if found in the aliases column.

    Parameters:
    - bigg_id (str): The BiGG ID to search for.
    - df (pd.DataFrame): The DataFrame containing 'aliases' and 'abbreviation' columns.

    Returns:
    - str: The corresponding KEGG ID from the abbreviation column, or 'NA' if not found.
    """
    for _, row in df.iterrows():
        aliases = row.get("aliases", [])
        if not isinstance(aliases, list):
            continue

        for alias in aliases:
            if alias.startswith("BiGG:"):
                # Split BiGG IDs and check for match
                ids = [i.strip() for i in alias.replace("BiGG:", "").split(";")]
                if bigg_id in ids:
                    return row.get("abbreviation", "NA")

    return "NA"


def fill_missing_kegg_ids_in_dict(mapping_dict, reactions_pandas):
    """
    Fill in missing KEGG IDs (None values) using a provided function that maps
    BiGG IDs to KEGG IDs.

    Parameters:
    - mapping_dict (dict): Dictionary with BiGG IDs as keys and KEGG IDs (or None) as values.
    - reactions_pandas (df): Dataframe from reactions.json file

    Returns:
    - dict: Updated dictionary with KEGG IDs filled where possible.
    """
    updated_dict = {}

    for bigg_id, kegg_id in mapping_dict.items():
        if kegg_id is not None:
            updated_dict[bigg_id] = kegg_id
        else:
            # Attempt to get KEGG ID using provided function
            guessed_kegg = bigg_to_kegg_id(bigg_id, reactions_pandas)
            updated_dict[bigg_id] = guessed_kegg

    return updated_dict


def get_kegg_pathways_from_reaction_ids(bigg_to_kegg_dict):
    """
    Given a dictionary mapping BiGG reaction IDs to KEGG reaction IDs,
    fetch the pathway information from KEGG for each reaction.

    Parameters:
    - bigg_to_kegg_dict (dict): Dictionary where keys are BiGG reaction IDs and values are KEGG reaction IDs.

    Returns:
    - pd.DataFrame: A DataFrame with columns ['bigg_reaction', 'kegg_reaction', 'pathway_ids', 'pathway_names'].
    """
    kegg = KEGG()
    reaction_to_pathways = {}

    for bigg_rxn, kegg_rxn in bigg_to_kegg_dict.items():
        if kegg_rxn:
            try:
                entry = kegg.get(kegg_rxn)
                parsed = kegg.parse(entry)
                pathways = parsed.get('PATHWAY', {})
                ids = list(pathways.keys())
                names = list(pathways.values())
                reaction_to_pathways[bigg_rxn] = {
                    'pathway_ids': ids,
                    'pathway_names': names
                }
            except Exception as e:
                print(f"Error retrieving {kegg_rxn}: {e}")
                reaction_to_pathways[bigg_rxn] = {
                    'pathway_ids': [],
                    'pathway_names': []
                }
        else:
            reaction_to_pathways[bigg_rxn] = {
                'pathway_ids': [],
                'pathway_names': []
            }

    # Convert to DataFrame
    df = pd.DataFrame([
        {
            'bigg_reaction': rxn,
            'kegg_reaction': bigg_to_kegg_dict[rxn],
            'pathway_ids': data['pathway_ids'],
            'pathway_names': data['pathway_names']
        }
        for rxn, data in reaction_to_pathways.items()
    ])

    return df


def subset_model_reactions_from_pathway_info(df, pathway_query):
    """
    Given a DataFrame with columns ['bigg_reaction', 'kegg_reaction', 'pathway_ids', 'pathway_names'],
    return all BiGG reaction IDs affiliated with a given KEGG pathway name or ID.

    Parameters:
    - df: pandas DataFrame from get_kegg_pathways_from_bigg_reactions.
    - pathway_query: str, KEGG pathway name (full or partial) or pathway ID (e.g., 'Glycolysis' or 'rn00010').

    Returns:
    - List of BiGG reaction IDs affiliated with the pathway.
    """
    query = pathway_query.strip().lower()

    matching_reactions = []

    for _, row in df.iterrows():
        ids = row.get('pathway_ids', [])
        names = row.get('pathway_names', [])

        # Normalize to lowercase for comparison
        ids_lower = [pid.lower() for pid in ids]
        names_lower = [pname.lower() for pname in names]

        # Check for match in either names or IDs
        if any(query in pid for pid in ids_lower) or any(query in pname for pname in names_lower):
            matching_reactions.append(row['bigg_reaction'])

    return sorted(set(matching_reactions))


def sort_reactions_in_pathways_by_reactions_in_model_order(full_list, *subsets):
    """
    Sorts and merges all subsets into a single list based on the order of elements in the full list.
    
    Parameters:
    - full_list (list): The reference list that defines the desired order.
    - *subsets (lists): One or more subset lists to be merged and ordered.
    
    Returns:
    - A single merged list of all subsets sorted by the order in full_list.
    
    Example Usage:
    
    glycolysis = ["PGI", "PFK", "FBA", "TPI", "GAPD", "PGK", "PGM", "ENO", "PYK"]
    ppp = ["G6PDH2r", "PGL", "GND", "RPE", "RPI", "TKT1", "TKT2", "TALA"]

    reactions_ordered = sort_reactions_in_pathways_by_reactions_in_model_order(ec_dingo_reactions, glycolysis, ppp)

    """
    index_map = {item: idx for idx, item in enumerate(full_list)}
    
    # Flatten all subsets into one list
    merged = [item for subset in subsets for item in subset]
    
    # Sort merged list by order in full_list
    sorted_merged = sorted(merged, key=lambda x: index_map.get(x, float('inf')))
    
    return sorted_merged


def subset_sampling_df_from_reaction_ids(samples, model_reactions, subset_reactions=[]):
    
    rxn_index_list = []
    
    for rxn in subset_reactions:
        rxn_index = model_reactions.index(rxn)
        rxn_index_list.append(rxn_index)
        
    subset_samples = samples[rxn_index_list]
    
    return subset_samples


def dictionary_bigg_id_to_pathway_names(**named_lists):
    bigg_id_to_pathway_dict = {}
    
    for name, list in named_lists.items():
        for item in list:
            if item not in bigg_id_to_pathway_dict:
                bigg_id_to_pathway_dict[item] = name
            elif bigg_id_to_pathway_dict[item] != name:
                bigg_id_to_pathway_dict[item] = "Multiple-Pathways"
    
    return bigg_id_to_pathway_dict


def dictionary_forward_reverse_bigg_id_to_pathway_names(bigg_to_pathway_dict, elements):
    bigg_id_forward_reverse_to_pathway_dict = {}
    
    for elem in elements:
        base_elem = elem
        if elem.endswith("_rev"):
            base_elem = elem[:-4]
        
        if base_elem in bigg_to_pathway_dict:
            bigg_id_forward_reverse_to_pathway_dict[elem] = bigg_to_pathway_dict[base_elem]
        else:
            bigg_id_forward_reverse_to_pathway_dict[elem] = "Unknown-Pathway"
    
    return bigg_id_forward_reverse_to_pathway_dict


def reaction_in_pathway_binary_matrix(mapping_dict):
    """
    Given a mapping dictionary, builds a binary matrix where rows are keys,
    columns are unique values, and the cell is 1 if the key maps to that value.
    """
    rows = list(mapping_dict.keys())
    cols = sorted(set(mapping_dict.values()))
    
    matrix = pd.DataFrame(0, index=rows, columns=cols)
    
    for key, value in mapping_dict.items():
        if value in matrix.columns:
            matrix.loc[key, value] = 1
            
    return matrix

   
def plot_reaction_in_pathway_heatmap(binary_df, font_size=12, fig_width=600, fig_height=400, title=""):
    """
    Plot a binary mapping matrix using Plotly with discrete color values (0 = grey, 1 = blue).
    
    Parameters:
    - binary_df: DataFrame with binary values (0 or 1)
    - font_size: Font size for axis labels and ticks
    - fig_width: Width of the figure in pixels
    - fig_height: Height of the figure in pixels
    """
    color_map = {
        0: "lightgrey",
        1: "steelblue"
    }

    z = binary_df.values
    x_labels = binary_df.columns
    y_labels = binary_df.index

    # Generate custom hover text per cell
    hover_text = [
        [f"Row: {y_labels[i]}<br>Col: {x_labels[j]}<br>Value: {z[i][j]}" for j in range(len(x_labels))]
        for i in range(len(y_labels))
    ]

    # Define discrete color scale
    colors = [[0, color_map[0]], [1, color_map[1]]]

    fig = go.Figure(
        data=go.Heatmap(
            z=z,
            x=x_labels,
            y=y_labels,
            hoverinfo="text",
            text=hover_text,
            colorscale=colors,
            showscale=False  # Hide the continuous colorbar
        )
    )

    # Add manual discrete legend
    for val, color in color_map.items():
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            marker=dict(size=10, color=color),
            name=str(val)
        ))

    fig.update_layout(
        title=title,
        xaxis=dict(tickangle=45, tickfont=dict(size=font_size)),
        yaxis=dict(tickfont=dict(size=font_size)),
        width=fig_width,
        height=fig_height,
        font=dict(size=font_size),
        legend_title="Value",
        margin=dict(t=40, b=40, l=40, r=40)
    )

    fig.show()