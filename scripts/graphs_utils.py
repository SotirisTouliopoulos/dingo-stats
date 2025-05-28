
import networkx as nx
import numpy as np
import plotly.graph_objects as go


def construct_graph(linear_correlation_matrix = None, non_linear_correlation_matrix = None, reactions=[], remove_unconnected_nodes=False, correction=True):
           
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


def plot_graph(G, pos):
    fig = go.Figure()

    # Edges
    for u, v, data in G.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        source = data.get('source', 'matrix1')
        weight = data.get('weight', 1.0)

        if weight is None or np.isnan(weight):
            weight = 0  # fallback to small value

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
        node_name = G.nodes[node].get('name', str(node))

        fig.add_trace(go.Scatter(
            x=[x], y=[y],
            mode='markers+text',
            marker=dict(size=10, color='black'),
            text=[node_name],
            textposition='top center',
            name=node_name,
            showlegend=False
        ))

    fig.update_layout(width=800, height=800, title='')
    fig.show()
    
    
def compare_betweenness_centralities(Graph_1, Graph_2):
    
    centrality_G1 = nx.betweenness_centrality(Graph_1, weight='weight', normalized=True)
    centrality_G2 = nx.betweenness_centrality(Graph_2, weight='weight', normalized=True)
    centrality_diff = {node: centrality_G2.get(node, 0) - centrality_G1.get(node, 0) for node in set(Graph_1.nodes()).union(Graph_2.nodes())}
    sorted_nodes = sorted(centrality_diff.items(), key=lambda x: x[1], reverse=True)

    return sorted_nodes
    #for node, change in sorted_nodes:
    #        print(f"Node {node}: ΔCentrality = {change:.4f}")