
import numpy as np
import plotly.express as px


def linear_correlations_in_reactions(samples, cobra_model, pearson_cutoff = 0.90):
    
    cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]

    # if provided sampling dataset has reactions as cols ==> transpose
    if samples.shape[1] == len(cobra_reactions_str):
        samples = samples.T

    correlation_matrix = np.corrcoef(samples)
    correlation_matrix[np.isnan(correlation_matrix)] = 0
    
    # find indices of correlation matrix where correlation does not occur
    no_corr_indices = np.argwhere((correlation_matrix < pearson_cutoff) & (correlation_matrix > -pearson_cutoff))
        
    filtered_correlation_matrix = correlation_matrix.copy()
    # replace values from the correlation matrix that do not overcome
    # the pearson cutoff with 0
    for i in range(0, no_corr_indices.shape[0]):        
        index1 = no_corr_indices[i][0]
        index2 = no_corr_indices[i][1]
        
        filtered_correlation_matrix[index1, index2] = 0
        
    return correlation_matrix, filtered_correlation_matrix


def plot_correlation_matrix(correlation_matrix, cobra_model):
    
    cobra_reactions_str = [str(reaction.id) for reaction in cobra_model.reactions]

    sns_colormap = [[0.0, '#3f7f93'],
                    [0.1, '#6397a7'],
                    [0.2, '#88b1bd'],
                    [0.3, '#acc9d2'],
                    [0.4, '#d1e2e7'],
                    [0.5, '#f2f2f2'],
                    [0.6, '#f6cdd0'],
                    [0.7, '#efa8ad'],
                    [0.8, '#e8848b'],
                    [0.9, '#e15e68'],
                    [1.0, '#da3b46']]
    
    fig = px.imshow(correlation_matrix, 
                    color_continuous_scale = sns_colormap,
                    x = cobra_reactions_str, y = cobra_reactions_str, origin="upper")
    
    fig.update_layout(
    xaxis=dict(tickfont=dict(size=5)),
    yaxis=dict(tickfont=dict(size=5)),
    width=900, height=900, plot_bgcolor="rgba(0,0,0,0)")
    
    fig.update_traces(xgap=1, ygap=1,   hoverongaps=False)
    
    fig.show()
