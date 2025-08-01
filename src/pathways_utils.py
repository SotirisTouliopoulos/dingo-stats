
from __future__ import annotations
import json
import pandas as pd
import plotly.graph_objects as go
from typing import Dict, Tuple, List, Optional
import numpy as np
from numpy.typing import NDArray
from Bio.KEGG import REST
from concurrent.futures import ThreadPoolExecutor, as_completed
from cobra import Model


def read_json_file(
    filepath: str
) -> Tuple[Dict, pd.DataFrame]:
    """
    Function that loads and returns the JSON file with KEGG information in JSON and pandas dataframe format
    
    Keyword arguments:
    filepath (str) -- path to the JSON file
    
    Tuple[dict, pd.DataFrame]
        reactions_json (dict) -- The raw content of the JSON file as a Python dictionary.
        reactions_pandas (pd.DataFrame) -- A pandas DataFrame constructed from the JSON content.
    """
    
    with open(filepath, 'r') as file:
        reactions_json = json.load(file)
        reactions_pandas = pd.DataFrame(reactions_json)

        return reactions_json, reactions_pandas
    

def map_model_to_kegg_reactions_dictionary(
    cobra_model: Model
) -> Dict[str, str]:
    """
    Function that creates a dictionary that will assign KEGG terms (values) to BiGG/SEED ids (keys) 
    only from model information (without searching in online databases or in external files)
    
    Keyword arguments:
    cobra_model (cobra.Model) -- cobra model object
    
    Returns:
    Dictionary (dict) -- maps each model's BiGG/SEED ID → KEGG ID 
    """
    
    initial_model_to_kegg_dictionary = {
        reac.id: reac.annotation.get("kegg.reaction") for reac in cobra_model.reactions}

    return initial_model_to_kegg_dictionary


def dictionary_reaction_id_to_kegg_id(
    reactions_pandas: pd.DataFrame
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Function that given the `reactions.json` file builds two dictionaries for fast lookup of KEGG reaction IDs from BiGG or SEED IDs.
    These dictionaries will be used as input in the `reaction_id_to_kegg_id` function.

    Keyword arguments:
    reactions_pandas (pd.DataFrame) -- DataFrame with columns including:
        'aliases' (list of strings): may contain entries like "BiGG:SUCDi" or "KEGG:R00010"
        'linked_reaction' (str): may contain SEED reaction IDs like "rxn12345;rxn67890"

    Returns:
    Tuple[dict, dict]:
        bigg_to_kegg (dict) -- maps each BiGG ID → KEGG ID
        seed_to_kegg (dict) -- maps each SEED ID → KEGG ID
    """

    # Helper to get the first KEGG alias from the 'aliases' list
    def first_kegg(aliases: Optional[List[str]]) -> Optional[str]:
        if not isinstance(aliases, list):
            return None
        for a in aliases:
            if isinstance(a, str) and a.startswith("KEGG:"):
                # Return just the KEGG ID part (e.g., "R00010")
                return a.split(":", 1)[1].strip()
        return None

    # Extract all BiGG IDs from the 'aliases' field
    def extract_bigg_ids(aliases: Optional[List[str]]) -> List[str]:
        if not isinstance(aliases, list):
            return []
        out: List[str] = []
        for a in aliases:
            if isinstance(a, str) and a.startswith("BiGG:"):
                # Split on semicolon, clean each ID
                out.extend([x.strip() for x in a.split(":", 1)[1].split(";") if x.strip()])
        return out

    # Extract all SEED IDs from the 'linked_reaction' string
    def extract_seed_ids(linked: Optional[str]) -> List[str]:
        if isinstance(linked, str):
            return [x.strip() for x in linked.split(";") if x.strip()]
        return []

    # Add temporary columns to the DataFrame for processing
    df = reactions_pandas.copy()

    # Column with KEGG ID (first KEGG tag found)
    df["kegg"]     = df["aliases"].apply(first_kegg)

    # Column with list of BiGG IDs extracted from aliases
    df["bigg_ids"] = df["aliases"].apply(extract_bigg_ids)

    # Column with list of SEED IDs extracted from linked_reaction
    df["seed_ids"] = df["linked_reaction"].apply(extract_seed_ids)

    # Initialize mapping dictionaries
    bigg_to_kegg: Dict[str, str] = {}
    seed_to_kegg: Dict[str, str] = {}

    # Fill bigg_to_kegg dictionary
    for bigg_ids, kegg in zip(df["bigg_ids"], df["kegg"]):
        if kegg is None:
            continue
        for bid in bigg_ids:
            # Only map the ID if not already added
            bigg_to_kegg.setdefault(bid, kegg)

    # Fill seed_to_kegg dictionary
    for seed_ids, kegg in zip(df["seed_ids"], df["kegg"]):
        if kegg is None:
            continue
        for sid in seed_ids:
            seed_to_kegg.setdefault(sid, kegg)

    return bigg_to_kegg, seed_to_kegg


def reaction_id_to_kegg_id(
    reaction_id: str,
    modeltype: str,
    bigg_to_kegg: Dict[str, str],
    seed_to_kegg: Dict[str, str],
) -> str:
    """
    FUnction that performs lookup to get KEGG ID for a given BiGG or SEED reaction ID.
    This function is used inside the `fill_missing_kegg_ids_in_initial_dictionary` function

    Keyword arguments:
    reaction_id (str) -- The BiGG or SEED reaction ID (e.g., "SUCDi" or "rxn12345")
    modeltype (str) -- Either "BiGG" or "SEED" (determines which dictionary to use)
    bigg_to_kegg (dict) -- Dictionary mapping BiGG IDs to KEGG IDs
    seed_to_kegg (dict) -- Dictionary mapping SEED IDs to KEGG IDs

    Returns:
    str -- The corresponding KEGG ID (e.g., "R00010"), or "NA" if not found
    """

    if modeltype == "BiGG":
        return bigg_to_kegg.get(reaction_id, "NA")
    elif modeltype == "SEED":
        return seed_to_kegg.get(reaction_id, "NA")
    else:
        raise ValueError(f"Unknown modeltype: {modeltype}")


def fill_missing_kegg_ids_in_initial_dictionary(
    initial_model_to_kegg_dictionary: dict[str, str], 
    bigg_to_kegg: Dict[str, str], 
    seed_to_kegg: Dict[str, str],
    modeltype: str = "BiGG"
) -> Dict:
    """
    Function that fills in missing KEGG IDs (NAs) in the initial mapping dictionary using a function that maps
    BiGG/SEED IDs to KEGG IDs.

    Keyword arguments:
    initial_model_to_kegg_dictionary (dict) -- Dictionary with reaction IDs as keys and KEGG IDs (or None) as values.
    bigg_to_kegg (dict) -- Dictionary mapping BiGG IDs to KEGG IDs
    seed_to_kegg (dict) -- Dictionary mapping SEED IDs to KEGG IDs
    modeltype (str) -- Either "BiGG" or "SEED" (determines which dictionary to use)
    
    Returns:
    final_model_to_kegg_dictionary (dict) -- Updated dictionary with KEGG IDs filled where possible.
    """
    
    final_model_to_kegg_dictionary = {}

    for reaction_id, kegg_id in initial_model_to_kegg_dictionary.items():
        if kegg_id is not None:
            final_model_to_kegg_dictionary[reaction_id] = kegg_id
        else:
            # Attempt to get KEGG ID using provided function
            mapped_kegg = reaction_id_to_kegg_id(reaction_id, modeltype, bigg_to_kegg, seed_to_kegg)
            final_model_to_kegg_dictionary[reaction_id] = mapped_kegg

    return final_model_to_kegg_dictionary


def fetch_kegg_pathway(
    kegg_rxn: str
) -> Tuple[str, List, List]:
    """
    Function that extracts associated KEGG pathways from a given KEGG id.

    Keyword arguments:
        kegg_rxn (str) -- KEGG reaction ID (e.g., 'R00010').
    
    Returns:
    Tuple[str, List, List]
        kegg_rxn (str) -- The KEGG reaction ID
        pathway_ids (list) -- The corresponding KEGG pathway IDs
        pathway_names (list) -- The corresponding KEGG pathway names
    """
    
    try:
        response = REST.kegg_get(kegg_rxn).read()
        lines = response.splitlines()
        pathways = {}
        in_pathway_section = False
        INDENT = " " * 12  # KEGG continuation lines are indented by 12 spaces

        for line in lines:
            if line.startswith("PATHWAY"):
                in_pathway_section = True
                # Everything after col 12: "rn00010  Glycolysis / Gluconeogenesis"
                tail = line[12:].strip()
                if tail:
                    pid, pname = tail.split(None, 1)  # split once on first whitespace
                    pathways[pid] = pname
            elif in_pathway_section and line.startswith(INDENT):
                tail = line[12:].strip()
                if tail:
                    pid, pname = tail.split(None, 1)
                    pathways[pid] = pname
            else:
                in_pathway_section = False

        return kegg_rxn, list(pathways.keys()), list(pathways.values())

    except Exception:
        return kegg_rxn, [], []


def get_kegg_pathways_from_reaction_ids(
    final_model_to_kegg_dictionary: Dict[str, str], 
    max_workers: int = 8
) -> pd.DataFrame:
    """
    Function that fetches KEGG pathway information for a set of model reactions using parallel requests.

    Keyword arguments:
    final_model_to_kegg_dictionary (dict) -- dictionary where keys are model reaction IDs and values are KEGG reaction IDs (e.g., 'R00010').
    max_workers (int) -- Number of threads for parallel downloading (default = 8).

    Returns:
    kegg_info_df (pd.DataFrame) -- A DataFrame with columns: ['model_reaction', 'kegg_reaction', 'pathway_ids', 'pathway_names']
    """

    # Filter out missing or empty KEGG IDs
    valid_kegg_rxns = {
        model_rxn: kegg_rxn.strip()
        for model_rxn, kegg_rxn in final_model_to_kegg_dictionary.items()
        if kegg_rxn and kegg_rxn.strip()
    }

    # Create a reverse mapping to track model_rxn by kegg_rxn
    reverse_map = {kegg_rxn: model_rxn for model_rxn, kegg_rxn in valid_kegg_rxns.items()}

    results = {}

    # Use threading to fetch KEGG entries in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(fetch_kegg_pathway, kegg_rxn): kegg_rxn
            for kegg_rxn in reverse_map
        }

        for future in as_completed(futures):
            kegg_rxn, pathway_ids, pathway_names = future.result()
            model_rxn = reverse_map[kegg_rxn]
            results[model_rxn] = {
                'kegg_reaction': kegg_rxn,
                'pathway_ids': pathway_ids,
                'pathway_names': pathway_names
            }

    # Fill in missing entries that failed or had no KEGG ID
    for model_rxn in final_model_to_kegg_dictionary:
        if model_rxn not in results:
            results[model_rxn] = {
                'kegg_reaction': final_model_to_kegg_dictionary[model_rxn],
                'pathway_ids': [],
                'pathway_names': []
            }

    # Convert results to DataFrame
    kegg_info_df = pd.DataFrame.from_dict(results, orient="index").reset_index()
    kegg_info_df.rename(columns={"index": "model_reaction"}, inplace=True)
    return kegg_info_df


def subset_model_reactions_from_pathway_info(
    kegg_info_df: pd.DataFrame, 
    pathway_query: str
) -> List:
    """
    Function that given a DataFrame with columns ['model_reaction', 'kegg_reaction', 'pathway_ids', 'pathway_names'],
    created wuth the `get_kegg_pathways_from_reaction_ids` function returns all reaction IDs affiliated 
    with a given KEGG pathway name or ID.

    Keyword arguments:
    kegg_info_df (pd.DataFrame) -- Output from `get_kegg_pathways_from_reaction_ids`, must contain 'pathway_ids' and 'pathway_names'.
    pathway_query (str) -- Exact KEGG pathway name or ID to match (e.g., 'Glycolysis / Gluconeogenesis' or 'rn00010').

    Returns:
    List[str] -- List of reaction IDs affiliated with the **exact** given pathway.
    """
    
    query = pathway_query.strip().lower()  # Normalize query for case-insensitive comparison

    matching_reactions = []

    for _, row in kegg_info_df.iterrows():
        # Extract and normalize pathway IDs and names
        ids = row.get('pathway_ids', [])
        names = row.get('pathway_names', [])

        ids_lower = [pid.lower() for pid in ids]
        names_lower = [pname.lower() for pname in names]

        # Match query exactly against any ID or name
        if query in ids_lower or query in names_lower:
            matching_reactions.append(row['model_reaction'])

    return sorted(set(matching_reactions))


def dictionary_reaction_id_to_pathway(
    **named_lists: List[str]
) -> Dict[str, str]:
    """
    Function that takes one or multiple lists containing reaction IDs (corresponding to KEGG pathways
    and creates a dictionary that maps the IDs to pathway names. If a reaction appears in more than 1 pathway,
    it is classified with the term "Multiple-Pathways"
    
    Keyword arguments:
    **named_lists: List[str] -- Named keyword arguments where each argument is a list of reaction IDs 
                                and the argument name represents the pathway name.
    
    Returns:
    reaction_id_to_pathway_dict (dict) -- dictionary mapping reaction id to pathway name
    """
    
    reaction_id_to_pathway_dict = {}
    
    for name, list in named_lists.items():
        for item in list:
            if item not in reaction_id_to_pathway_dict:
                reaction_id_to_pathway_dict[item] = name
            elif reaction_id_to_pathway_dict[item] != name:
                reaction_id_to_pathway_dict[item] = "Multiple-Pathways"
    
    return reaction_id_to_pathway_dict


def reaction_in_pathway_binary_matrix(
    reaction_id_to_pathway_dict: Dict
) -> pd.DataFrame:
    """
    Function that given a mapping dictionary, builds a binary matrix where rows (reactions) are keys,
    columns (pathways) are unique values, and the cell is 1 if the key maps to that value.
    
    Keyword arguments:
    mapping_dict (Dict) -- dictionary mapping reaction id to pathway name
    
    Returns:
    binary_df (pd.DataFrame) -- DataFrame with binary values (0 or 1) matching reactions to pathways
    """

    rows = list(reaction_id_to_pathway_dict.keys())
    cols = sorted(set(reaction_id_to_pathway_dict.values()))

    binary_df = pd.DataFrame(0, index=rows, columns=cols)

    for key, value in reaction_id_to_pathway_dict.items():
        if value in binary_df.columns:
            binary_df.loc[key, value] = 1

    return binary_df


def plot_reaction_in_pathway_heatmap(
    binary_df: pd.DataFrame,
    font_size: int = 12,
    fig_width: int = 600,
    fig_height: int = 400,
    title: str = ""
):
    """
    Function that plots a binary mapping matrix created from the `reaction_in_pathway_binary_matrix` function.
    
    Keyword arguments:
    binary_df (pd.DataFrame) -- DataFrame with binary values (0 or 1)
    font_size (int) -- Font size for axis labels and ticks
    fig_width (int) -- Width of the figure in pixels
    fig_height (int) -- Height of the figure in pixels
    title (str) -- Title of the plot
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


def sort_reactions_by_model_order(
    full_list: List, 
    *subsets: List
) -> List:
    """
    Function that flattens the lists provided in the `subsets` argument (corresponding to reactions from different pathways) 
    in a single list and then orders the element of the new list based on the order of the reaction in the initial model

    Keyword arguments:
    full_list (List) -- The reference list that defines the desired order.
    *subsets (List) -- One or more subset lists to be merged and ordered.
    
    Returns:
    sorted_merged -- a single merged list of all subsets sorted by the order in full_list.
    
    Example Usage:
    Glycolysis = ["PGI", "PFK", "FBA", "TPI", "GAPD", "PGK", "PGM", "ENO", "PYK"]
    PPP = ["G6PDH2r", "PGL", "GND", "RPE", "RPI", "TKT1", "TKT2", "TALA"]
    reactions_ordered = sort_reactions_in_pathways_by_reactions_in_model_order(ec_cobra_reactions_str, Glycolysis, PPP)
    """
    
    index_map = {item: idx for idx, item in enumerate(full_list)}
    
    # Flatten all subsets into one list
    merged = [item for subset in subsets for item in subset]
    
    # Sort merged list by order in full_list
    sorted_merged = sorted(merged, key=lambda x: index_map.get(x, float('inf')))
    
    return sorted_merged


def subset_sampling_array_from_reaction_ids(
    samples: NDArray[np.float64], 
    model_reactions: List, 
    subset_reactions: List =[]
) -> NDArray[np.float64]:
    """
    Function that takes a sampling array with reactions as rows and samples as columns and subsets it
    to include only reactions of interest
    
    Keyword arguments:
    samples (Numpy 2D array) -- A sampling 2D array with reactions as rows and samples as columns
    model_reactions (List) -- A list containing the model's reactions
    subset_reactions (List) -- A list containing reactions of interest to subset the sampling array
    
    Returns:
    subset_samples (NDArray[np.float64]) -- subset of the sampling dataframe containing only reactions of interest
    """
    
    rxn_index_list = []
    
    for rxn in subset_reactions:
        rxn_index = model_reactions.index(rxn)
        rxn_index_list.append(rxn_index)
        
    subset_samples = samples[rxn_index_list]
    
    return subset_samples


def dictionary_map_reverse_reaction_id_to_pathway(
    reaction_id_to_pathway_dict: Dict, 
    for_rev_reactions: List
) -> Dict[str, str]:
    """
    Function that is used when we split bidirectional reactions to separate forward and reverse reactions.
    It maps the reverse reaction to the corresponding pathway (the one that the forward reactions maps to)
    It enriches the dictionary created from the `dictionary_reaction_id_to_pathway` function
    
    Keyword arguments:
    reaction_id_to_pathway_dict (Dict) -- Dict mapping reaction IDs to pathway names
    for_rev_reactions (List) -- List of the splitted reactions
    
    Returns:
    reaction_id_forward_reverse_to_pathway_dict (Dict) -- Dictionary containing reaction-pathway mapping information 
                                                          for forward and reverse reactions separately
    """
    
    reaction_id_forward_reverse_to_pathway_dict = {}
    
    for elem in for_rev_reactions:
        base_elem = elem
        if elem.endswith("_rev"):
            base_elem = elem[:-4]
        
        if base_elem in reaction_id_to_pathway_dict:
            reaction_id_forward_reverse_to_pathway_dict[elem] = reaction_id_to_pathway_dict[base_elem]
        else:
            reaction_id_forward_reverse_to_pathway_dict[elem] = "Other-Pathways"
    
    return reaction_id_forward_reverse_to_pathway_dict