
import networkx as nx
import numpy as np
import plotly.graph_objects as go
from networkx.algorithms.community import modularity, greedy_modularity_communities
import plotly.express as px
import matplotlib.cm as cm
import matplotlib.colors as mcolors


def construct_graph(linear_correlation_matrix=None, 
                    non_linear_correlation_matrix=None, 
                    reactions=[], 
                    remove_unconnected_nodes=False, 
                    correction=True, 
                    group_map=None):
    """
    Creates a graph from a linear correlation matrix or from both a linear and a non-linear correlation matrix.
    Adds node attributes for group information if group_map is provided.

    Parameters:
    - linear_correlation_matrix: numpy 2D array
    - non_linear_correlation_matrix: numpy 2D array (optional)
    - reactions: list of reaction names (ordered like matrix indices)
    - remove_unconnected_nodes: if True, removes isolated nodes
    - correction: if True, use absolute values of correlations
    - group_map: dictionary mapping reaction names to group names (optional)
      
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
    - NetworkX graph
    - Layout positions for plotting
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


def plot_graph(G, pos):
    """Plot graph with nodes colored by group and a legend."""
    
    fig = go.Figure()

    # Edges
    for u, v, data in G.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        source = data.get('source', 'matrix1')
        weight = data.get('weight', 1.0)
        
        #edge_color = 'blue' if data['weight'] > 0 else 'red'
        
        if source == 'matrix1':
            if weight >= 0:
                color = 'blue'
                dash = 'solid'
            else:
                color = 'red'
                dash = 'solid'
        elif source == 'matrix2':
            if weight >= 0:
                color = 'blue'
                dash = 'dash'
            else:
                color = 'red'
                dash = 'dot'
        else:
            color = 'gray'
            dash = 'solid'

        
        fig.add_trace(go.Scatter(
            x=[x0, x1], y=[y0, y1], mode='lines',
            line=dict(width=abs(data['weight']) * 1, color=color, dash=dash),
            hoverinfo='none',
            showlegend=False
        ))

    # Node group collection
    group_nodes = {}
    for node in G.nodes():
        group = G.nodes[node].get('group', 'Unknown')
        group_nodes.setdefault(group, []).append(node)
    
    # Color map
    color_map = {group: f'rgb{tuple((np.array(cm.tab10(i % 10)[:3]) * 255).astype(int))}' 
                 for i, group in enumerate(sorted(group_nodes.keys()))}
    
    # Nodes by group
    for group, nodes in group_nodes.items():
        x_vals, y_vals, labels = [], [], []
        for node in nodes:
            x, y = pos[node]
            x_vals.append(x)
            y_vals.append(y)
            labels.append(node)
        
        fig.add_trace(go.Scatter(
            x=x_vals, y=y_vals, mode='markers+text',
            marker=dict(size=10, color=color_map[group]),
            text=labels,
            textposition='top center',
            name=group,
            showlegend=True
        ))

    fig.update_layout(
        width=800, height=800,
        legend=dict(title='Groups')
    )
    fig.show()
    
    
def compare_network_modularities(Graph_1, Graph_2, tol=1e-3):
    """
    Function that compares the modularity between 2 different graphs
    """
    
    for u, v, d in Graph_1.edges(data=True):
        d['weight'] = abs(d['weight'])
        
    for u, v, d in Graph_2.edges(data=True):
        d['weight'] = abs(d['weight'])
    
    communities_graph_1 = greedy_modularity_communities(Graph_1, weight="weight")
    mod_score_graph_1 = modularity(Graph_1, communities_graph_1, weight="weight")
    
    communities_graph_2 = greedy_modularity_communities(Graph_2, weight="weight")
    mod_score_graph_2 = modularity(Graph_2, communities_graph_2, weight="weight")
    
    print(mod_score_graph_1, mod_score_graph_2)
    
    mod_diff = mod_score_graph_2 - mod_score_graph_1
    
    if mod_diff > tol:
        print("Increased Modularity: Graph_2")
    elif mod_diff < tol:
        print("Increased Modularity: Graph_1")
    else:
        print("Graphs with similar Modularity")
    
    return mod_diff


def compare_betweenness_centralities(Graph_1, Graph_2):
    """
    Function that compares the betweenness centralities across all nodes between 2 different graphs
    
    Correlations must be converted, as this function uses higher values as distances
    """
    
    for u, v, d in Graph_1.edges(data=True):
        d['distance'] = 1 - abs(d['weight'])
        
    for u, v, d in Graph_2.edges(data=True):
        d['distance'] = 1 - abs(d['weight'])
    
    centrality_G1 = nx.betweenness_centrality(Graph_1, weight='distance', normalized=True)
    centrality_G2 = nx.betweenness_centrality(Graph_2, weight='distance', normalized=True)
    
    all_nodes = set(Graph_1.nodes()).union(Graph_2.nodes())

    centrality_diff = {
        node: (centrality_G2.get(node, 0) - centrality_G1.get(node, 0))
        for node in all_nodes
    }
    
    sorted_nodes = sorted(centrality_diff.items(), key=lambda x: x[1], reverse=True)

    return sorted_nodes


def compare_closeness_centralities(Graph_1, Graph_2):
    """
    Function that compares the closeness centralities across all nodes between 2 different graphs
    
    Correlations must be converted, as this function uses higher values as distances
    """
    
    for u, v, d in Graph_1.edges(data=True):
        d['distance'] = 1 - abs(d['weight'])
        
    for u, v, d in Graph_2.edges(data=True):
        d['distance'] = 1 - abs(d['weight'])
        
    centrality_G1 = nx.closeness_centrality(Graph_1, distance='distance')
    centrality_G2 = nx.closeness_centrality(Graph_2, distance='distance')
    
    all_nodes = set(Graph_1.nodes()).union(Graph_2.nodes())
    centrality_diff = {
        node: centrality_G2.get(node, 0) - centrality_G1.get(node, 0)
        for node in all_nodes
    }
    
    sorted_nodes = sorted(centrality_diff.items(), key=lambda x: x[1], reverse=True)
    return sorted_nodes


def weighted_node_centrality(Graph_1, Graph_2):
    """
    Function that creates the difference between node centralities (sum of values of weights) between 2 different graphs
    """
    
    for u, v, d in Graph_1.edges(data=True):
        d['weight'] = abs(d['weight'])
        
    for u, v, d in Graph_2.edges(data=True):
        d['weight'] = abs(d['weight'])
        
    all_nodes = set(Graph_1.nodes()).union(Graph_2.nodes())

    centrality_diff = {}
    for node in all_nodes:
        deg1 = Graph_1.degree(node, weight='weight') if node in Graph_1 else 0
        deg2 = Graph_2.degree(node, weight='weight') if node in Graph_2 else 0
        centrality_diff[node] = deg2 - deg1

    sorted_nodes = sorted(centrality_diff.items(), key=lambda x: x[1], reverse=True)
    return sorted_nodes