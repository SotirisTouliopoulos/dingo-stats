
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import math
from scipy.spatial.distance import jensenshannon


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


def plot_correlation_matrix(correlation_matrix, reactions=[]):
    
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
                    x = reactions, y = reactions, origin="upper")
    
    fig.update_layout(
    xaxis=dict(tickfont=dict(size=5)),
    yaxis=dict(tickfont=dict(size=5)),
    width=900, height=900, plot_bgcolor="rgba(0,0,0,0)")
    
    fig.update_traces(xgap=1, ygap=1,   hoverongaps=False)
    
    fig.show()


def compute_copula(flux1, flux2, n):
    """A Python function to estimate the copula between two fluxes

    Keyword arguments:
    flux1: A vector that contains the measurements of the first reaxtion flux
    flux2: A vector that contains the measurements of the second reaxtion flux
    n: The number of cells
    """

    N = flux1.size
    copula = np.zeros([n,n], dtype=float)

    I1 = np.argsort(flux1)
    I2 = np.argsort(flux2)

    grouped_flux1 = np.zeros(N)
    grouped_flux2 = np.zeros(N)

    for j in range(n):
        rng = range((j*math.floor(N/n)),((j+1)*math.floor(N/n)))
        grouped_flux1[I1[rng]] = j
        grouped_flux2[I2[rng]] = j
    
    for i in range(n):
        for j in range(n):
            copula[i,j] = sum((grouped_flux1==i) *( grouped_flux2==j))
    
    copula = copula / N
    return copula


def plot_copula(data_flux1, data_flux2, n = 5, width = 900 , height = 600, export_format = "svg"):
    """A Python function to plot the copula between two fluxes

    Keyword arguments:
    data_flux1: A list that contains: (i) the vector of the measurements of the first reaction,
                                      (ii) the name of the first reaction
    data_flux2: A list that contains: (i) the vector of the measurements of the second reaction,
                                      (ii) the name of the second reaction
    n: The number of cells
    """

    flux1 = data_flux1[0]
    flux2 = data_flux2[0]
    copula = compute_copula(flux1, flux2, n)

    fig = go.Figure(
            data   = [go.Surface(z=copula)],
            layout = go.Layout(
                height = height, 
                width  = width,
            )
        )


    fig.update_layout(
            title = 'Copula between '+ data_flux1[1] + ' and ' + data_flux2[1],
            scene = dict(
                    xaxis_title= data_flux1[1],
                    yaxis_title= data_flux2[1],
                    zaxis_title="prob, mass"
                ),
            margin=dict(r=30, b=30, l=30, t=50))

    fig.layout.template = None
    
    fig.show()

    camera = dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=1.25, y=1.25, z=1.25)
    )

    fig.update_layout(scene_camera=camera)
    fig.to_image(format = export_format, engine="kaleido")
    
    

def split_forward_reverse(steady_states, reactions=[]):
    extended_reactions = []
    extended_steady_states = []
    for index, r in enumerate(steady_states):
        if (r < 0).any():
            s = -(r)
            extended_steady_states.append(s)
            extended_reactions.append("_".join([reactions[index], "rev"]))
                
        extended_steady_states.append(r)
        extended_reactions.append(reactions[index])
        
    extended_steady_states = np.array(extended_steady_states)
    
    return extended_steady_states, extended_reactions
        

def correlated_reactions(steady_states, split_bidirectional_reactions=False, include_non_linear=False, reactions=[], 
                         pearson_cutoff = 0.30, indicator_cutoff = 0.5, jensenshannon_cutoff = 0.1, std_cutoff = 1e-6,
                         cells = 10, cop_coeff = 0.3, 
                         lower_triangle = True, verbose = False):
    """A Python function to calculate the correlation matrix from a steady states array

    Keyword arguments:
    steady_states -- A numpy array of the generated steady states fluxes
    split_bidirectional_reactions -- A boolean variable that if True, splits flux from reversible reactions into forward and reverse flux
    include_non_linear -- A boolean variable that if True, takes into account and calculates non-linear correlations
    reactions -- A list with the reactions IDs (must be in accordance with the rows of the steady states)
    pearson_cutoff -- A cutoff to filter (remove) linear correlations based on the pearson coefficient
    jensenshannon_cutoff -- A cutoff to filter (remove) non-linear correlations based on the Jensen-Shannon metric
    std_cutoff -- A cutoff to avoid computing the copula between 2 fluxes with almost fixed values
    indicator_cutoff -- A cutoff to classify non-linear correlations as positive, negative or non-significant
    cells -- Number of cells to compute the copula
    cop_coeff -- A value that narrows or widens the width of the copula's diagonal
    lower_triangle -- A boolean variable that if True returns only the lower triangular matrix
    verbose -- A boolean variable that if True additional information is printed as an output.
    
    Returns:
    linear_correlation_matrix -- A correlations matrix filtered based on the given cutoffs that includes only linear correlations
    mixed_correlation_matrix -- A correlations matrix filtered based on the given cutoffs that includes both linear and non-linear correlations
    non_linear_correlation_matrix -- A correlations matrix filtered based on the given cutoffs that includes only non-linear correlations
    correlations_dictionary -- A dictionary containing unique reaction pairs and their corresponding correlation values
    """
    
    if indicator_cutoff < 1:
        raise Exception("Indicator cutoff must be at least equal to 1")
    
    if split_bidirectional_reactions == False:
        steady_states = np.absolute(steady_states)
    else:
        steady_states, reactions = split_forward_reverse(steady_states, reactions)
        
    
    if cop_coeff > 0.4 or cop_coeff < 0.2:
        raise Exception("Input value to cop_coeff parameter must be between 0.2 and 0.4")
    
    # calculate coefficients to access red and blue copula mass
    cop_coeff_1 = cop_coeff
    cop_coeff_2 = 1 - cop_coeff
    cop_coeff_3 = 1 + cop_coeff
    
    # compute correlation matrix
    corr_matrix = np.corrcoef(steady_states, rowvar=True)
    # replace not assigned values with 0
    corr_matrix[np.isnan(corr_matrix)] = 0
    # keep only the lower triangle (with the diagonal) to reduce computational time
    corr_matrix[np.triu_indices(corr_matrix.shape[0], 1)] = np.nan

    # create a copy of correlation matrix to replace/filter values
    linear_correlation_matrix = corr_matrix.copy()    
    
    # find indices of correlation matrix where correlation does not occur
    no_corr_indices = np.argwhere((linear_correlation_matrix < pearson_cutoff) & (linear_correlation_matrix > -pearson_cutoff))
    # find indices of correlation matrix where correlation does occur
    corr_indices = np.argwhere((linear_correlation_matrix > pearson_cutoff) | (linear_correlation_matrix < -pearson_cutoff))
        
    # replace values from the correlation matrix that do not overcome the pearson cutoff with 0
    for i in range(0, no_corr_indices.shape[0]):       
        index1 = no_corr_indices[i][0]
        index2 = no_corr_indices[i][1]
        
        if index1 == index2:
            continue
        else:
            linear_correlation_matrix[index1, index2] = 0
            linear_correlation_matrix[index2, index1] = 0  
        
        
    correlations_dictionary = {}
    
    for i in range(0, corr_indices.shape[0]):       
        index1 = corr_indices[i][0]
        index2 = corr_indices[i][1]
        
        if index1 == index2:
            continue
        else:
            reaction1 = reactions[index1]
            reaction2 = reactions[index2]
            
            pearson = linear_correlation_matrix[index1, index2]
            
            if pearson > 0:
                correlations_dictionary[reaction1 + "~" + reaction2] = {'pearson': pearson,
                                                                  'jensenshannon': 0,
                                                                  'indicator': 0,
                                                                  'classification': "positive"}
            elif pearson < 0:
                correlations_dictionary[reaction1 + "~" + reaction2] = {'pearson': pearson,
                                                                  'jensenshannon': 0,
                                                                  'indicator': 0,
                                                                  'classification': "negative"}
                

    # if user does not want to calculate non-linear correlations
    if include_non_linear == False:
        linear_correlation_matrix[np.triu_indices(linear_correlation_matrix.shape[0], 1)] = np.nan
        np.fill_diagonal(linear_correlation_matrix, 1)

        if lower_triangle == True:
            
            return linear_correlation_matrix, correlations_dictionary
        
        else:
            # fill the upper triangle and return a square correlation matrix
            linear_correlation_matrix = np.tril(linear_correlation_matrix) + np.tril(linear_correlation_matrix, -1).T
            return linear_correlation_matrix, correlations_dictionary
    
    else:
        # A correlation matrix that will store both linear and and non-linear correlations
        mixed_correlation_matrix = linear_correlation_matrix.copy()
        
        # compute copula for each set of non-correlated reactions
        for i in range(0, no_corr_indices.shape[0]):
            
            index1 = no_corr_indices[i][0]
            index2 = no_corr_indices[i][1]
            
            if index1 == index2:
                continue
            
            reaction1 = reactions[index1]
            reaction2 = reactions[index2]
                        
            flux1 = steady_states[index1]
            flux2 = steady_states[index2]
            
            if (np.std(flux1) < std_cutoff) or (np.std(flux2) < std_cutoff):
                continue
                
            copula = compute_copula(flux1, flux2, cells)
            copula_flat = copula.flatten()

            rows, cols = copula.shape
        
            red_mass = 0
            blue_mass = 0
            indicator = 0
                            
            for row in range(rows):
                for col in range(cols):
                    # values in the diagonal
                    if ((row-col >= -cop_coeff_1*rows) & (row-col <= cop_coeff_1*rows)):        
                        # values near the top left and bottom right corner
                        if ((row+col < cop_coeff_2*rows) | (row+col > cop_coeff_3*rows)):
                            red_mass = red_mass + copula[row][col]
                    else:
                        # values near the top right and bottom left corner
                        if ((row+col >= cop_coeff_2*rows-1) & (row+col <= cop_coeff_3*rows-1)):
                            blue_mass = blue_mass + copula[row][col]

            indicator = (red_mass+1e-9) / (blue_mass+1e-9)
            
                            
            # define a uniform copula to compare with previously computed copula
            reference_copula = np.full( (cells, cells), (1 / (cells*cells)) )
            reference_copula_flat = reference_copula.flatten()
            
            # calculate distance (scaled: 0-1) between copula from reaction pairs and reference copula
            js = jensenshannon(copula_flat, reference_copula_flat)                
                
            
            # classify specific pair of reactions as positive or negative correlated based on indicator cutoff 
            if (indicator > indicator_cutoff) and (js > jensenshannon_cutoff):                
                pearson = mixed_correlation_matrix[index1, index2]
                mixed_correlation_matrix[index1, index2] = abs(js)

                correlations_dictionary[reaction1 + "~" + reaction2] = {'pearson': pearson,
                                                                  'jensenshannon': js,
                                                                  'indicator': indicator,
                                                                  'classification': "positive"}

            elif (indicator < 1/indicator_cutoff) and (js > jensenshannon_cutoff):
                pearson = mixed_correlation_matrix[index1, index2]

                mixed_correlation_matrix[index1, index2] = -abs(js)
                
                correlations_dictionary[reaction1 + "~" + reaction2] = {'pearson': pearson,
                                                                  'jensenshannon': js,
                                                                  'indicator': indicator,
                                                                  'classification': "negative"}
                
            else:
                pearson = mixed_correlation_matrix[index1, index2]

                correlations_dictionary[reaction1 + "~" + reaction2] = {'pearson': pearson,
                                                                  'jensenshannon': js,
                                                                  'indicator': indicator,
                                                                  'classification': "no_correlation"}
                
            if verbose == True:
                print("Completed the process of",i+1,"from",no_corr_indices.shape[0],"copulas")
        
        
        # Fill upper triangle and diagonal with NAs
        mixed_correlation_matrix[np.triu_indices(mixed_correlation_matrix.shape[0], 1)] = np.nan
        np.fill_diagonal(mixed_correlation_matrix, 1)
        
        # Create a  matrix containing only non-linear correlations
        mask = (mixed_correlation_matrix != 0) & (linear_correlation_matrix == 0)
        non_linear_correlation_matrix = np.zeros_like(mixed_correlation_matrix)
        non_linear_correlation_matrix[mask] = mixed_correlation_matrix[mask]


        if lower_triangle == True:
            return linear_correlation_matrix, non_linear_correlation_matrix, mixed_correlation_matrix, correlations_dictionary

        else:
            # fill the upper triangle and return a square correlation matrix
            mixed_correlation_matrix = np.tril(mixed_correlation_matrix) + np.tril(mixed_correlation_matrix, -1).T
            return linear_correlation_matrix, non_linear_correlation_matrix, mixed_correlation_matrix, correlations_dictionary
