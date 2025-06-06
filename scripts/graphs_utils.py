
import networkx as nx
import numpy as np
import plotly.graph_objects as go
from networkx.algorithms.community import modularity, greedy_modularity_communities
import plotly.express as px
import matplotlib.cm as cm
import matplotlib.colors as mcolors


def construct_graph(linear_correlation_matrix = None, non_linear_correlation_matrix = None, reactions=[], remove_unconnected_nodes=False, correction=True):
    """
    Function used to create a graph either from a linear correlation matrix, or from both a linear and non-linear correlation matrix 
    """      
    
    graph_matrix_1 = linear_correlation_matrix.copy()
    np.fill_diagonal(graph_matrix_1, 0)
    graph_matrix_1 = np.nan_to_num(graph_matrix_1, nan=0.0)
    
    if correction == True:
        graph_matrix_1 = np.absolute(graph_matrix_1)
        
    G = nx.from_numpy_array(graph_matrix_1)
    
    for u, v, d in G.edges(data=True):
        d['source'] = 'matrix1'
        
    for u, v, d in G.edges(data=True):
        G.add_edge(u, v, weight=d['weight'], source='matrix1')
                
    
    if non_linear_correlation_matrix is not None:
        G_combined = G.copy()

        graph_matrix_2 = non_linear_correlation_matrix.copy()
        np.fill_diagonal(graph_matrix_2, 0)    
        graph_matrix_2 = np.nan_to_num(graph_matrix_2, nan=0.0)
        
        if correction == True:
            graph_matrix_2 = np.absolute(graph_matrix_2)
            
        G2 = nx.from_numpy_array(graph_matrix_2)
        
        for u, v, d in G2.edges(data=True):
            d['source'] = 'matrix2'
            
        for u, v, d in G2.edges(data=True):
            if G_combined.has_edge(u, v):
                continue
            G_combined.add_edge(u, v, weight=d['weight'], source='matrix2')
            
        G_combined = nx.relabel_nodes(G_combined, lambda x: reactions[x] if reactions else x)
        
        #pos = nx.spring_layout(G_combined)
        unconnected_nodes = list(nx.isolates(G_combined))
        if remove_unconnected_nodes == True:
            G_combined.remove_nodes_from(unconnected_nodes)
        
        pos = nx.spring_layout(G_combined)
    
        return G_combined, pos
    
    
    G = nx.relabel_nodes(G, lambda x: reactions[x] if reactions else x)
    
    unconnected_nodes = list(nx.isolates(G))
    if remove_unconnected_nodes == True:
        G.remove_nodes_from(unconnected_nodes)
    
    pos = nx.spring_layout(G)
    
    return G, pos


def plot_graph(G, pos, clusters, showlegend=False):
    """
    Function that given a graph created from the function above, plots its layout in an interaxtive plot.
    
    clusters is a nested list showing the structure of the clusters created from the "cluster_corr_reactions" of dingo.
    When showlegend is True, nodes appear on the right of the panel and user can turn off nodes of no interest
    """
    
    fig = go.Figure()
    
    # Define color palette
    def get_continuous_colors(n):
        cmap = cm.get_cmap("viridis", n)
        return [mcolors.to_hex(cmap(i)) for i in range(n)]

    cluster_colors = get_continuous_colors(len(clusters))

    # Assign cluster index to each node
    node_cluster_map = {}
    for idx, cluster in enumerate(clusters):
        for node in cluster:
            node_cluster_map[node] = idx

    # Edges
    for u, v, data in G.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        source = data.get('source', 'matrix1')
        weight = data.get('weight', 1.0)

        if weight is None or np.isnan(weight):
            weight = 0

        # Fixed color by source
        if source == 'matrix1':
            color = 'blue' if weight >= 0 else 'lightblue'
        elif source == 'matrix2':
            color = 'red' if weight >= 0 else 'orange'
        else:
            color = 'gray'

        # Edge width scales with correlation strength
        line_width = abs(weight) * 1  # adjust multiplier as needed

        fig.add_trace(go.Scatter(
            x=[x0, x1], y=[y0, y1],
            mode='lines',
            line=dict(width=line_width, color=color),
            hoverinfo='text',
            text=f'{u} ↔ {v}<br>Weight: {weight:.2f}<br>Source: {source}',
            showlegend=False
        ))

    # Nodes
    for node in G.nodes():
        x, y = pos[node]
        cluster_idx = node_cluster_map.get(node, -1)
        node_color = cluster_colors[cluster_idx % len(cluster_colors)] if cluster_idx >= 0 else 'black'
        node_name = G.nodes[node].get('name', str(node))

        fig.add_trace(go.Scatter(
            x=[x], y=[y],
            mode='markers+text',
            marker=dict(size=10, color=node_color),
            text=[node_name],
            textposition='top center',
            name=node_name,
            showlegend=showlegend
        ))

    fig.update_layout(width=1000, height=1000, title='')
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