
import networkx as nx
import numpy as np
import plotly.graph_objects as go
from networkx.algorithms.community import modularity, greedy_modularity_communities
import plotly.express as px
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from scipy.spatial import ConvexHull
from itertools import combinations
from networkx.algorithms.clique import find_cliques
import cobra
from typing import Dict, Tuple, List, Any, Optional
from cobra import Model, Reaction
from numpy.typing import NDArray


def construct_graph(
    linear_correlation_matrix: NDArray[np.float64] = None, 
    non_linear_correlation_matrix: NDArray[np.float64] = None, 
    reactions: List = [], 
    remove_unconnected_nodes: bool = False, 
    correction: bool = True, 
    group_map: Dict = None
) -> Tuple[nx.Graph, Dict[Any, Tuple[float, float]]]:
    """
    Function that creates a networkx graph from a linear correlation matrix or from both a linear correlation 
    and a non-linear copula dependencies matrix. In this graph reactions are nodes and correlation values are edges.
    Users can also provide reaction-pathway mapping information in the `group_map` dictionary parameter 
    and this will be saved in the graph for potential vizualization.

    Keyword arguments:
    linear_correlation_matrix (NDArray[np.float64]) -- numpy 2D array corresponding to the linear correlation matrix
    non_linear_correlation_matrix (NDArray[np.float64]) numpy 2D array (optional) corresonding to the non-linear copula dependencies matrix
    reactions (List) -- list of reaction names (ordered like matrix indices)
    remove_unconnected_nodes (bool) -- if True, removes isolated nodes
    correction (bool) -- if True, use absolute values of correlations
    group_map (Dict) -- dictionary mapping reaction names to group names (pathways)
      
        group_map = {
                "PGI" : "Glycolysis",
                "PGI_rev" : "Glycolysis", 
                "PFK" : "Glycolysis", 
                ...,
                "G6PDH2r" : "Pentose-Phosphate", 
                "PGL" : "Pentose-Phosphate", 
                "GND" : "Pentose-Phosphate", 
                }

    Returns:
    Tuple[nx.Graph, Dict[Any, Tuple[float, float]]
    G/G_combined (nx.Graph) -- The created NetworkX graph
    pos (Dict[Any, Tuple[float, float]]) -- dictionary of node positions, usually from a layout function
    """
    
    graph_matrix_1 = linear_correlation_matrix.copy()
    np.fill_diagonal(graph_matrix_1, 0)
    graph_matrix_1 = np.nan_to_num(graph_matrix_1, nan=0.0)
    
    if correction:
        graph_matrix_1 = np.abs(graph_matrix_1)
        
    G = nx.from_numpy_array(graph_matrix_1)
    
    for u, v, d in G.edges(data=True):
        d['source'] = 'matrix1'
        d['weight'] = graph_matrix_1[u, v]
        
    if non_linear_correlation_matrix is not None:
        G_combined = G.copy()

        graph_matrix_2 = non_linear_correlation_matrix.copy()
        np.fill_diagonal(graph_matrix_2, 0)
        graph_matrix_2 = np.nan_to_num(graph_matrix_2, nan=0.0)
        
        if correction:
            graph_matrix_2 = np.abs(graph_matrix_2)
        
        G2 = nx.from_numpy_array(graph_matrix_2)
        
        for u, v, d in G2.edges(data=True):
            d['source'] = 'matrix2'
            d['weight'] = graph_matrix_2[u, v]
            
        for u, v, d in G2.edges(data=True):
            if not G_combined.has_edge(u, v):
                G_combined.add_edge(u, v, weight=d['weight'], source='matrix2')
        
        # Relabel nodes
        G_combined = nx.relabel_nodes(G_combined, lambda x: reactions[x] if reactions else x)

        # Add group info if available
        if group_map:
            for node in G_combined.nodes():
                G_combined.nodes[node]['group'] = group_map.get(node, "Unknown")

        if remove_unconnected_nodes:
            unconnected_nodes = list(nx.isolates(G_combined))
            G_combined.remove_nodes_from(unconnected_nodes)
        
        pos = nx.spring_layout(G_combined)
        return G_combined, pos
    
    # Only linear matrix used
    G = nx.relabel_nodes(G, lambda x: reactions[x] if reactions else x)

    if group_map:
        for node in G.nodes():
            G.nodes[node]['group'] = group_map.get(node, "Unknown")

    if remove_unconnected_nodes:
        unconnected_nodes = list(nx.isolates(G))
        G.remove_nodes_from(unconnected_nodes)

    pos = nx.spring_layout(G)
    return G, pos


def draw_positive_clique_shadows(
    G: nx.Graph,
    pos: Dict[Any, Tuple[float, float]],
    fig: go.Figure,
    positive_cliques: List[List[Any]],
    min_clique_size: int = 3
):
    """
    Function that adds translucent shaded areas around positive cliques to a Plotly figure. This function highlights 
    regions in a network where nodes form strongly positively connected cliques by drawing convex hull 
    polygons around them.

    Keyword arguments:
    G (nx.Graph) -- NetworkX graph containing nodes and weighted edges
    pos (Dict) -- dictionary of node positions, usually from a layout function
    fig (go.Figure) -- Plotly figure object to which shaded regions will be added
    positive_cliques (List[List]) -- list of positive cliques, where each clique is a list of node names
    min_clique_size (int) -- minimum number of nodes in a clique for it to be visualized (default = 3)
    """
    
    fillcolor = 'rgba(144, 238, 144, 0.5)'  # light green

    for clique in positive_cliques:
        if len(clique) < min_clique_size:
            continue

        # Confirm all pairwise edges in the clique are positive
        is_positive = True
        for i, u in enumerate(clique):
            for v in clique[i+1:]:
                if not G.has_edge(u, v) or G[u][v]['weight'] <= 0:
                    is_positive = False
                    break
            if not is_positive:
                break

        if not is_positive:
            continue

        pts = np.array([pos[n] for n in clique])
        if len(pts) >= 3:
            try:
                hull = ConvexHull(pts)
                hull_pts = pts[hull.vertices]
                x_hull = list(hull_pts[:, 0]) + [hull_pts[0, 0]]
                y_hull = list(hull_pts[:, 1]) + [hull_pts[0, 1]]

                fig.add_trace(go.Scatter(
                    x=x_hull,
                    y=y_hull,
                    fill='toself',
                    mode='lines',
                    line=dict(color='rgba(0,0,0,0)', width=0),
                    fillcolor=fillcolor,
                    hoverinfo='skip',
                    showlegend=False,
                    name='Positive Clique Shadow'
                ))
            except Exception as e:
                print(f"ConvexHull error for positive clique shadow: {e}")


def draw_negative_clique_shadows(
    G: nx.Graph,
    pos: Dict[Any, Tuple[float, float]],
    fig: go.Figure,
    positive_cliques: List[List[Any]],
    min_clique_size: int = 3
):
    """
    Function that adds translucent shaded areas around regions where two positive cliques are connected by at least one 
    negative edge. This function visually highlights antagonistic relationships between strongly connected 
    positive subgraphs by drawing convex hulls around their union.

    Keyword arguments:
    G (nx.Graph) -- NetworkX graph containing nodes and weighted edges
    pos (Dict[Any, Tuple[float, float]]) -- dictionary of node positions, usually from a layout function
    fig (go.Figure) -- Plotly figure object to which shaded regions will be added
    positive_cliques (List[List[Any]]) -- list of positive cliques, where each clique is a list of node names
    min_clique_size (int) -- minimum number of nodes in a clique for it to be considered (default = 3)
    """
    
    fillcolor = 'rgba(255, 0, 0, 0.15)'  # light red

    for clique1, clique2 in combinations(positive_cliques, 2):
        if len(clique1) < min_clique_size or len(clique2) < min_clique_size:
            continue

        # Check for any negative edge between the cliques
        has_negative_edge = any(
            G.has_edge(u, v) and G[u][v]['weight'] < 0
            for u in clique1 for v in clique2
        )

        if not has_negative_edge:
            continue

        combined_nodes = clique1 + clique2
        pts = np.array([pos[n] for n in combined_nodes])

        if len(pts) >= 3:
            try:
                hull = ConvexHull(pts)
                hull_pts = pts[hull.vertices]
                x_hull = list(hull_pts[:, 0]) + [hull_pts[0, 0]]
                y_hull = list(hull_pts[:, 1]) + [hull_pts[0, 1]]

                fig.add_trace(go.Scatter(
                    x=x_hull,
                    y=y_hull,
                    fill='toself',
                    mode='lines',
                    line=dict(color='rgba(0,0,0,0)', width=0),
                    fillcolor=fillcolor,
                    hoverinfo='skip',
                    showlegend=False,
                    name='Negative Clique Shadow'
                ))
            except Exception as e:
                print(f"ConvexHull error for negative clique shadow: {e}")

              
def plot_graph(
    G: nx.Graph,
    pos: Dict[Any, Tuple[float, float]],
    centralities: Dict[str, Dict[Any, float]],
    remove_clique_edges: bool = False,
    include_matrix2_in_cliques: bool = True,
    min_clique_size: int = 5,
    shadow_edges: Optional[str] = None
):
    """
    Function that plots a correlation-based networkx graph with multiple visual features including node annotations,
    edge styles, and clique-based shadowing.

    Keyword arguments:
    G (nx.Graph) -- NetworkX graph with nodes and weighted edges; edges should have 'weight' and 'source' attributes
    pos (Dict[Any, Tuple[float, float]]) -- dictionary of node positions, usually from a layout function
    centralities (Dict[str, Dict[Any, float]]) -- dictionary containing precomputed centrality metrics
        Example:
            {
                "degree": {node1: val1, ...},
                "betweenness": {node1: val2, ...},
                "clustering": {node1: val3, ...}
            }
    remove_clique_edges (bool) -- if True, remove edges within positive cliques and between cliques with negative connections
    include_matrix2_in_cliques (bool) -- if False, only use matrix1 (linear edges) for positive clique detection
    min_clique_size (int) -- minimum size for cliques to be considered (default = 5)
    shadow_edges (Optional[str]) -- visualization mode for clique shadows:
        - None: no shadows
        - "positive": show shadows around positive cliques
        - "negative": show shadows between cliques with negative edges
        - "mixed": show both types of shadows
    """

    fig = go.Figure()

    # Subgraph for positive clique detection
    if include_matrix2_in_cliques:
        G_pos = G.edge_subgraph([(u, v) for u, v, d in G.edges(data=True) if d['weight'] > 0]).copy()
    else:
        G_pos = G.edge_subgraph([
            (u, v) for u, v, d in G.edges(data=True)
            if d['weight'] > 0 and d.get('source') == 'matrix1'
        ]).copy()

    positive_cliques = [clique for clique in find_cliques(G_pos) if len(clique) >= min_clique_size]

    # Track positive clique nodes
    positive_clique_nodes = set()
    for clique in positive_cliques:
        positive_clique_nodes.update(clique)

    # Annotate all nodes with neighbors and clique membership
    for node in G.nodes():
        neighbors = list(G.neighbors(node))
        pos_neighbors = [v for v in neighbors if G[node][v]['weight'] > 0]
        neg_neighbors = [v for v in neighbors if G[node][v]['weight'] < 0]

        G.nodes[node]['positive_neighbors'] = sorted(pos_neighbors)
        G.nodes[node]['negative_neighbors'] = sorted(neg_neighbors)
        G.nodes[node]['in_clique'] = node in positive_clique_nodes

    # Prepare graph copy for plotting
    G_to_plot = G.copy()

    if remove_clique_edges:
        # Remove internal positive clique edges
        edges_to_remove = set()
        for clique in positive_cliques:
            for i in range(len(clique)):
                for j in range(i + 1, len(clique)):
                    u, v = clique[i], clique[j]
                    if G_to_plot.has_edge(u, v) and G_to_plot[u][v]['weight'] > 0:
                        edges_to_remove.add((u, v))
        G_to_plot.remove_edges_from(edges_to_remove)

        # Remove negative edges between cliques
        for clique1, clique2 in combinations(positive_cliques, 2):
            neg_edges = [
                (u, v) for u in clique1 for v in clique2
                if G_to_plot.has_edge(u, v) and G_to_plot[u][v]['weight'] < 0
            ]
            G_to_plot.remove_edges_from(neg_edges)

        if shadow_edges == None:
            pass
        elif shadow_edges == "positive":
            draw_positive_clique_shadows(G_pos, pos, fig, positive_cliques, min_clique_size)
        elif shadow_edges == "negative":
            draw_negative_clique_shadows(G, pos, fig, positive_cliques, min_clique_size)
        elif shadow_edges == "mixed":
            draw_positive_clique_shadows(G_pos, pos, fig, positive_cliques, min_clique_size)
            draw_negative_clique_shadows(G, pos, fig, positive_cliques, min_clique_size)
            

    # Group edges by style for efficient plotting
    edge_groups = {}
    for u, v, data in G_to_plot.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        weight = data.get('weight', 1.0)
        source = data.get('source', 'matrix1')

        if weight >= 0:
            if source == 'matrix1':
                style_key = 'Linear Positive'
                color = 'blue'
                dash = 'solid'
            elif source == 'matrix2':
                style_key = 'Non-linear Positive'
                color = 'blue'
                dash = 'dash'
            else:
                style_key = 'Unknown Positive'
                color = 'gray'
                dash = 'solid'
        else:
            if source == 'matrix1':
                style_key = 'Linear Negative'
                color = 'red'
                dash = 'solid'
            elif source == 'matrix2':
                style_key = 'Non-linear Negative'
                color = 'red'
                dash = 'dot'
            else:
                style_key = 'Unknown Negative'
                color = 'gray'
                dash = 'solid'

        if style_key not in edge_groups:
            edge_groups[style_key] = {'x': [], 'y': [], 'color': color, 'dash': dash}

        edge_groups[style_key]['x'] += [x0, x1, None]
        edge_groups[style_key]['y'] += [y0, y1, None]

    for style_key, edge_data in edge_groups.items():
        fig.add_trace(go.Scatter(
            x=edge_data['x'],
            y=edge_data['y'],
            mode='lines',
            line=dict(width=2, color=edge_data['color'], dash=edge_data['dash']),
            hoverinfo='none',
            name=style_key,
            showlegend=True
        ))

    # Group nodes by 'group' attribute and assign colors
    group_nodes = {}
    for node in G_to_plot.nodes():
        group = G_to_plot.nodes[node].get('group', 'Unknown')
        group_nodes.setdefault(group, []).append(node)

    color_map = {
    group: f"rgb{tuple(int(v) for v in (np.array(cm.tab10(i % 10)[:3]) * 255))}"
    for i, group in enumerate(sorted(group_nodes.keys()))
    }

    for group, nodes in group_nodes.items():
        x_vals, y_vals, labels, sizes = [], [], [], []
        for node in nodes:
            x, y = pos[node]
            x_vals.append(x)
            y_vals.append(y)

            pos_neighbors = G.nodes[node]['positive_neighbors']
            neg_neighbors = G.nodes[node]['negative_neighbors']
            in_clique = G.nodes[node]['in_clique']
             
            degree = centralities.get("degree").get(node, 0)
            betweenness = centralities.get("betweenness").get(node, 0)
            clustering = centralities.get("clustering").get(node, 0)

            # Truncate long neighbor lists
            def truncate_list(lst, max_len=10):
                return lst[:max_len] + ['...'] if len(lst) > max_len else lst

            pos_neighbors_trunc = truncate_list(pos_neighbors)
            neg_neighbors_trunc = truncate_list(neg_neighbors)

            label = (
                f"Node: {node}<br>"
                f"Pos: {', '.join(pos_neighbors_trunc) if pos_neighbors_trunc else 'None'}<br>"
                f"Neg: {', '.join(neg_neighbors_trunc) if neg_neighbors_trunc else 'None'}<br>"
                f"In clique: {in_clique}<br>"
                f"Degree centrality: {degree:.2f}<br>"
                f"Betweenness centrality: {betweenness:.2f}<br>"
                f"Clustering coefficient: {clustering:.2f}"
            )
            labels.append(label)
            sizes.append(14 if in_clique else 10)

        fig.add_trace(go.Scatter(
            x=x_vals, y=y_vals, mode='markers+text',
            marker=dict(size=sizes, color=color_map[group], line=dict(width=1, color='black')),
            text=[str(n) for n in nodes],
            textposition='top center',
            hovertext=labels,
            hoverinfo='text',
            name=group,
            showlegend=True,
        ))

    fig.update_layout(
        width=900,
        height=900,
        legend=dict(title='Edge Type & Node Group')
    )

    fig.show()


def compute_nodes_centrality_metrics(
    G: nx.Graph
) -> Dict[str, Dict[Any, float]]:
    """
    Function that computes centrality measures on a correlation-based graph.

    Keyword arguments:
    G (nx.Graph) -- NetworkX graph with edge weights representing correlation values in the range [-1, 1]

    Returns:
    Dict[str, Dict[Any, float]] -- A dictionary containing centrality metrics for each node:
        - 'degree': Weighted degree centrality normalized by number of nodes
        - 'betweenness': Betweenness centrality using distance = 1 - abs(weight)
        - 'clustering': Clustering coefficient using abs(weight) as edge weight
    """
    
    G_abs = G.copy()
    
    # Convert correlations to distances
    for u, v, d in G_abs.edges(data=True):
        abs_w = abs(d.get('weight', 0))
        
        d['abs_weight'] = float(abs_w)
        d['distance'] = float(1 - abs_w)
        
    # Compute weighted degree (sum of absolute edge weights)
    weighted_degree = {
        n: sum(G_abs[u][v]['abs_weight'] for u, v in G_abs.edges(n))
        for n in G_abs.nodes()
    }
    
    # Normalize weighted degree
    number_of_nodes = len(G_abs.nodes())
    weighted_degree_normalized = {
        k: v / (number_of_nodes - 1) if number_of_nodes > 1 else 0.0
        for k, v in weighted_degree.items()
    }
    
    # Betweenness centrality
    betweenness = nx.betweenness_centrality(G_abs, weight='distance')
    betweenness = {k: float(v) for k, v in betweenness.items()}

    # Clustering coefficient
    clustering = nx.clustering(G_abs, weight='abs_weight')
    clustering = {k: float(v) for k, v in clustering.items()}
    

    centrality_dict = {
        'degree': weighted_degree_normalized,
        'betweenness': betweenness,
        'clustering': clustering
    }
    
    return centrality_dict


def compare_essential_to_network_central_reactions(
    cobra_model: cobra.Model,
    centrality_dict: Dict[str, float],
    threshold: float = 0.999
) -> Tuple[int, int, List[str]]:    
    """
    Function that counts how many central network reactions (reactions with high centrality metrics)
    belong to the group of essential reactions 

    Keyword arguments:
    cobra_model (Model) -- COBRA metabolic model object
    centrality_dict (Dict[str, float]) -- dictionary mapping reaction IDs to centrality scores
    threshold (float) -- threshold to determine essentiality; a reaction is essential if its removal
                         causes the objective function to drop below `threshold * optimal value` (default = 0.999)

    Returns:
    Tuple[int, int, List[str]]
        Number of reactions that are both central and essential
        Total number of essential reactions
        Sorted list of reaction IDs that are both central and essential
    """
    
    def filter_reversible_reactions(
        top_reactions: List[Tuple[str, float]]
    ) -> List[Tuple[str, float]]:
        """
        FUnction that handle nodes from reversible reactions by keeping only one representative 
        (the one with the highest centrality score) for each forward/reverse pair. This is to avoid
        counting the same reversible reaction as central 2 times (the forward and the reverse)

        Keyword arguments:
        top_reactions (List[Tuple[str, float]]) -- list of (reaction_id, centrality_score) tuples, 
                                                   sorted in descending order of centrality

        Returns:
        List[Tuple[str, float]] -- filtered list with one entry per reaction (ignoring reverse duplicates)
        """

        seen = set()
        filtered = []

        for rxn_id, score in top_reactions:
            base_id = rxn_id.replace("_rev", "")

            if base_id not in seen:
                seen.add(base_id)
                filtered.append((base_id, score))

        return filtered
    
    solution = cobra_model.optimize()
    obj_sol = solution.objective_value
    obj_sol_threshold = obj_sol * threshold
    
    essential_reactions = cobra.flux_analysis.variability.find_essential_reactions(cobra_model, threshold=obj_sol_threshold)
    essential_reactions = [reaction.id for reaction in essential_reactions]
    
    essential_reactions_count = len(essential_reactions)

    # Sort reactions by centrality descending
    top_reactions = sorted(centrality_dict.items(), key=lambda x: x[1], reverse=True)
    # Keep the first node id (highest centrality value) from a forward and reverse pair
    top_reactions = filter_reversible_reactions(top_reactions)
    top_ids = set([r[0] for r in top_reactions[:essential_reactions_count]])

    # Convert reaction list to set
    input_set = set(essential_reactions)

    # Find intersection
    central_and_essential_reactions = input_set & top_ids
    central_and_essential_reactions_count = len(central_and_essential_reactions)
    
    print(top_reactions)

    return central_and_essential_reactions_count, essential_reactions_count, sorted(central_and_essential_reactions)


def compare_network_modularities(
    Graph_1: nx.Graph, 
    Graph_2: nx.Graph
) -> float:
    """
    Function that compares the modularity scores of two graphs using the greedy modularity community detection method
    and returns the difference substracting the modularity score of graph 2 (second argument) from graph 1 (first argument).

    Keyword arguments:
    Graph_1 (nx.Graph) -- first NetworkX graph
    Graph_2 (nx.Graph) -- second NetworkX graph
    
    Returns:
    float -- difference in modularity (modularity of Graph_2 minus modularity of Graph_1)
    """
    
    for u, v, d in Graph_1.edges(data=True):
        d['weight'] = abs(d['weight'])
        
    for u, v, d in Graph_2.edges(data=True):
        d['weight'] = abs(d['weight'])
    
    communities_graph_1 = greedy_modularity_communities(Graph_1, weight="weight")
    mod_score_graph_1 = modularity(Graph_1, communities_graph_1, weight="weight")
    
    communities_graph_2 = greedy_modularity_communities(Graph_2, weight="weight")
    mod_score_graph_2 = modularity(Graph_2, communities_graph_2, weight="weight")
        
    mod_diff = mod_score_graph_2 - mod_score_graph_1
    return mod_diff


def compare_node_centralities(
    centrality_dict_1: Dict[str, float], 
    centrality_dict_2: Dict[str, float]
) -> List[Tuple[str, float]]:
    """
    Function that compares node centralities between two centrality dictionaries for shared nodes and returns 
    the difference by substracting the corresponding metric score of graph 2 (second argument) from graph 1 (first argument).

    Keyword arguments:
    centrality_dict_1 (Dict[str, float]) -- Dictionary mapping node IDs to centrality values (first graph)
    centrality_dict_2 (Dict[str, float]) -- Dictionary mapping node IDs to centrality values (second graph)

    Returns:
    Tuple -- A list of tuples (node, centrality_difference) sorted by difference descending,
    where centrality_difference = centrality_dict_2[node] - centrality_dict_1[node]
    """
    
    all_nodes = set(centrality_dict_1.keys()) & set(centrality_dict_2.keys())
    
    centrality_diff = {
        node: (centrality_dict_2.get(node, 0) - centrality_dict_1.get(node, 0))
        for node in all_nodes
    }
    
    sorted_nodes = sorted(centrality_diff.items(), key=lambda x: x[1], reverse=True)
    return sorted_nodes