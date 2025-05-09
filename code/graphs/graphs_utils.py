
import networkx as nx
import numpy as np
import plotly.graph_objects as go


def construct_graph(correlation_matrix, cobra_model, remove_unconnected_nodes=False, correction=True):
       
    cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]
    
    graph_matrix = correlation_matrix.copy()
    np.fill_diagonal(graph_matrix, 0)
    
    if correction == True:
        graph_matrix = abs(graph_matrix)
        
    G = nx.from_numpy_array(graph_matrix)
    G = nx.relabel_nodes(G, lambda x: cobra_reactions_str[x])
    
    pos = nx.spring_layout(G)
    unconnected_nodes = list(nx.isolates(G))

    if remove_unconnected_nodes == True:
        G.remove_nodes_from(unconnected_nodes)
    
    return G, pos


def plot_graph(G, pos):
    
    fig = go.Figure()

    for u, v, data in G.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        
        edge_color = 'blue' if data['weight'] > 0 else 'red'
        
        fig.add_trace(go.Scatter(x=[x0, x1], y=[y0, y1], mode='lines', 
                                    line=dict(width=abs(data['weight']) * 1, 
                                    color=edge_color), hoverinfo='none',
                                    showlegend=False))

    for node in G.nodes():
        x, y = pos[node]
        node_name = G.nodes[node].get('name', f'Node {node}')
  
        fig.add_trace(go.Scatter(x=[x], y=[y], mode='markers', 
                                    marker=dict(size=10),
                                    text=[node_name],
                                    textposition='top center',
                                    name = node_name,
                                    showlegend=False))
        
    fig.update_layout(width=800, height=800)
    fig.show()