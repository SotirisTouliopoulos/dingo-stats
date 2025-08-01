
import numpy as np
import plotly.graph_objects as go
from scipy.cluster import hierarchy
from matplotlib import cm
from collections import defaultdict
from typing import Dict, Tuple, List
from cobra import Model, Reaction
from numpy.typing import NDArray


def clustering_of_correlation_matrix(
    correlation_matrix: NDArray[np.float64], 
    reactions: List, 
    linkage: str = "ward",
    t: float = 4.0, 
    correction: bool = True
) -> Tuple[NDArray[np.float64], NDArray[np.float64], List[List]]:
    """
    Function for hierarchical clustering of the correlation matrix

    Keyword arguments:
    correlation_matrix (NDArray[np.float64]) -- A numpy 2D array of a correlation matrix
    reactions (List) -- A list with the corresponding reactions (must be in accordance with the correlation matrix)
    linkage (str) -- linkage defines the type of linkage. Available linkage types are: single, average, complete, ward.
    t (float) -- A threshold that defines a threshold that cuts the dendrogram at a specific height and produces clusters
    correction (bool) -- A boolean variable that if True converts the values of the the correlation matrix to absolute values.
                         
    Returns:
    Tuple[NDArray[np.float64], NDArray[np.float64], List[List]]
        dissimilarity_matrix (NDArray[np.float64]) -- The dissimilarity matrix calculated from the given correlation matrix
        labels (NDArray[np.float64]) -- Integer labels corresponding to clusters
        clusters (List[List]) -- Nested list containing reaction IDs grouped based on their cluster labels
    """
    
    def clusters_list(
        reactions: List, 
        labels: NDArray[np.float64],
        ) -> List[List]:
        """
        Function to return a nested list with grouped reactions based on clustering
        
        Keyword arguments:
        reactions (List) -- A list with the corresponding reactions (must be in accordance with the correlation matrix)
        labels (NDArray[np.float64]) -- Integer labels corresponding to clusters

        Returns: 
        clusters (List[List]) -- Nested list containing reaction IDs grouped based on their cluster labels
        """
        
        clusters = []
        unique_labels = np.unique(labels)
        for label in unique_labels:
            cluster = []
            label_where = np.where(labels == label)[0]
            for where in label_where:
                cluster.append(reactions[where])
            clusters.append(cluster)
        return clusters
    
    if correction == True:        
        dissimilarity_matrix = 1 - abs(correlation_matrix)
    else:
        dissimilarity_matrix = 1 - correlation_matrix
            
    Z = hierarchy.linkage(dissimilarity_matrix, linkage)
    labels = hierarchy.fcluster(Z, t, criterion='distance')
    
    clusters = clusters_list(reactions, labels)
    return dissimilarity_matrix, labels, clusters


def plot_dendrogram(
    dissimilarity_matrix: NDArray[np.float64],
    reactions: List,
    group_map: dict,
    linkage: str = 'ward',
    t: float = 4.0,
    label_fontsize: int = 10,
    height: int = 600,
    width: int = 1000,
    show_labels: bool = True,
    title: str = "" 
):
    """
    Function that plots the dendrogram from hierarchical clustering
    
    Keyword arguments:
    dissimilarity_matrix (NDArray[np.float64]) -- The dissimilarity matrix calculated from the `clustering_of_correlation_matrix` function
    reactions (List) -- List with reactions IDs
    group_map (Dict) -- Dictionary mapping reactions to pathways
    linkage (str) -- linkage defines the type of linkage. Available linkage types are: single, average, complete, ward.
    t (float) -- A threshold that defines a threshold that cuts the dendrogram at a specific height
    label_fontsize (int) -- Size for the plotted labels (reaction IDs)
    height (int) -- Defines the height of the plot
    width (int) -- Defines the width of the plot
    show_labels (bool) -- Boolean variable that if True labels are shown in the x-axis
    title (str) -- Title for the plot
    """
    
    # Default fallback
    if group_map is None:
        group_map = {}

    # If group_map is empty, assign all to 'Other'
    assigned_groups = [group_map.get(r, 'Other') for r in reactions]
    unique_groups = sorted(set(assigned_groups))
    color_map = {
        g: f'rgba{tuple((np.array(cm.tab10(i % 10))[:3] * 255).astype(int)) + (1.0,)}'
        for i, g in enumerate(unique_groups)
    }
    
    # Hierarchical clustering
    Z = hierarchy.linkage(dissimilarity_matrix, method=linkage)
    dendro = hierarchy.dendrogram(Z, labels=reactions, color_threshold=t, no_plot=True)

    # Extract coordinates of dendrogram lines
    icoord = np.array(dendro['icoord'])
    dcoord = np.array(dendro['dcoord'])
    labels = dendro['ivl']
    # leaves = dendro['leaves']

    # Group-color mapping
    if not group_map:
        group_map = {r: 'Other' for r in reactions}

    unique_groups = sorted(set(group_map.values()))
    
    color_map = {
    g: f'rgba({r},{g_},{b},1.0)'
    for i, g in enumerate(unique_groups)
    for r, g_, b in [tuple(int(x * 255) for x in cm.tab10(i % 10)[:3])]
    }
    
    # label_colors = {label: color_map[group_map.get(label, unique_groups[0])] for label in labels}

    # Plot dendrogram branches
    fig = go.Figure()

    for xs, ys in zip(icoord, dcoord):
        fig.add_trace(go.Scatter(
            x=xs, y=ys,
            mode='lines',
            line=dict(color='black', width=1),
            hoverinfo='none',
            showlegend=False
        ))

    # Plot colored leaf points
    leaf_labels = dendro['ivl']
    num_leaves = len(leaf_labels)
    leaf_xs = list(range(5, 10 * num_leaves + 5, 10))  # Matches default dendrogram spacing
    leaf_ys = [0] * num_leaves
    # leaf_colors = [label_colors[label] for label in leaf_labels]

    # Build per-group traces for toggle-able legend
    group_to_points = defaultdict(list)
    for x, y, label in zip(leaf_xs, leaf_ys, leaf_labels):
        group = group_map.get(label, 'Other')
        group_to_points[group].append((x, y, label))

    for group, points in group_to_points.items():
        xs, ys, labels_in_group = zip(*points)
        
        fig.add_trace(go.Scatter(
        x=xs,
        y=ys,
        mode='markers+text' if show_labels else 'markers',
        text=labels_in_group if show_labels else None,
        textposition='bottom center',
        textfont=dict(size=label_fontsize),
        marker=dict(color=color_map[group], size=8),
        name=group,
        legendgroup=group,
        showlegend=True
    ))

    # Draw threshold line
    fig.add_shape(
        type='line',
        x0=min(min(icoord.flatten()), 0),
        x1=max(max(icoord.flatten()), 0),
        y0=t,
        y1=t,
        line=dict(dash='dash', color='gray', width=1)
    )

    fig.update_layout(
        title=title,
        xaxis=dict(showticklabels=False),
        yaxis=dict(title='Distance'),
        width=width,
        height=height,
        legend_title="Groups"
    )

    fig.show()