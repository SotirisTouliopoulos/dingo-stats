
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import math
from scipy.spatial.distance import jensenshannon


def correlated_reactions(steady_states, reactions=[],
                         split_bidirectional_reactions=False, include_non_linear=False, boolean_sharing_metabolites_matrix=None, 
                         pearson_cutoff = 0.30, indicator_cutoff = 1.0, jensenshannon_cutoff = 0.1, std_cutoff = 1e-6,
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
    cop_coeff -- A value that narrows or widens the width of the copula's diagonal (use lower values to capture extreme tail dependences)
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
    # fill diagonal with 1
    np.fill_diagonal(corr_matrix, 1)
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
                classification =  "positive"
            elif pearson < 0:
                classification = "negative"

            correlations_dictionary[reaction1 + "~" + reaction2] = {'pearson': pearson,
                                                                'jensenshannon': 0,
                                                                'indicator': 0,
                                                                'classification': classification}
    
    # remove correlations from reactions not sharing metabolites if user provides the boolean matrix
    if boolean_sharing_metabolites_matrix is not None:
        linear_correlation_matrix[~boolean_sharing_metabolites_matrix] = 0
    
        
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
            dependence, indicator = copula_tail_dependence(copula, cop_coeff_1, cop_coeff_2, cop_coeff_3, indicator_cutoff)
            copula_flat = copula.flatten()
                            
            # define a uniform copula to compare with previously computed copula
            reference_copula = np.full( (cells, cells), (1 / (cells*cells)) )
            reference_copula_flat = reference_copula.flatten()
            
            # calculate distance (scaled: 0-1) between copula from reaction pairs and reference copula
            js = jensenshannon(copula_flat, reference_copula_flat)                
                
            
            # classify the dependenc ebetween a specific pair of reactions based on indicator metric
            if (js > jensenshannon_cutoff):
                if dependence.startswith("positive"):
                    js = abs(js)               
                    mixed_correlation_matrix[index1, index2] = js
                
                elif dependence.startswith("negative"):
                    js = -abs(js)
                    mixed_correlation_matrix[index1, index2] = js
                
                else:
                    js = 0
                    mixed_correlation_matrix[index1, index2] = js

                
                
                correlations_dictionary[reaction1 + "~" + reaction2] = {'pearson': 0,
                                                                    'jensenshannon': js,
                                                                    'indicator': indicator,
                                                                    'classification': dependence}
                
            
            if verbose == True:
                print("Completed the process of",i+1,"from",no_corr_indices.shape[0],"copulas")
        
        
        # Fill upper triangle and diagonal with NAs
        mixed_correlation_matrix[np.triu_indices(mixed_correlation_matrix.shape[0], 1)] = np.nan
        np.fill_diagonal(mixed_correlation_matrix, 1)
        
        # Create a  matrix containing only non-linear correlations
        mask = (mixed_correlation_matrix != 0) & (linear_correlation_matrix == 0)
        non_linear_correlation_matrix = np.zeros_like(mixed_correlation_matrix)
        non_linear_correlation_matrix[mask] = mixed_correlation_matrix[mask]
        
        # remove correlations from reactions not sharing metabolites if user provides the boolean matrix
        if boolean_sharing_metabolites_matrix is not None:
            non_linear_correlation_matrix[~boolean_sharing_metabolites_matrix] = 0
            linear_correlation_matrix[~boolean_sharing_metabolites_matrix] = 0
            mixed_correlation_matrix[~boolean_sharing_metabolites_matrix] = 0
            

        if lower_triangle == True:
            return linear_correlation_matrix, non_linear_correlation_matrix, mixed_correlation_matrix, correlations_dictionary

        else:
            # fill the upper triangle and return a square correlation matrix
            non_linear_correlation_matrix = np.tril(non_linear_correlation_matrix) + np.tril(non_linear_correlation_matrix, -1).T
            linear_correlation_matrix = np.tril(linear_correlation_matrix) + np.tril(linear_correlation_matrix, -1).T
            mixed_correlation_matrix = np.tril(mixed_correlation_matrix) + np.tril(mixed_correlation_matrix, -1).T
            
            return linear_correlation_matrix, non_linear_correlation_matrix, mixed_correlation_matrix, correlations_dictionary


def plot_correlation_matrix(correlation_matrix, reactions=[]):
    # Ensure data is clipped between -1 and +1
    correlation_matrix = np.clip(correlation_matrix, -1, 1)

    sns_colormap = [
        [0.0, '#d73027'],   # red for -1
        [0.5, '#f7f7f7'],   # white for 0
        [1.0, '#4575b4']    # blue for +1
    ]
    
    fig = px.imshow(
        correlation_matrix,
        color_continuous_scale=sns_colormap,
        zmin=-1,  # Fixed scale
        zmax=1,
        x=reactions,
        y=reactions,
        origin="upper"
    )
    
    fig.update_layout(
        xaxis=dict(tickfont=dict(size=5)),
        yaxis=dict(tickfont=dict(size=5)),
        width=900,
        height=900,
        plot_bgcolor="rgba(0,0,0,0)"
    )
    
    fig.update_traces(xgap=1, ygap=1, hoverongaps=False)
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


def copula_tail_dependence(copula, cop_coeff_1, cop_coeff_2, cop_coeff_3, indicator_cutoff):

    rows, cols = copula.shape
        
    red_mass = 0
    blue_mass = 0
    
    top_left = 0
    bottom_right = 0
    top_right = 0
    bottom_left = 0
    
    indicator = 0
                    
    for row in range(rows):
        for col in range(cols):
            # values in the diagonal
            if ((row-col >= -cop_coeff_1*rows) & (row-col <= cop_coeff_1*rows)): 
                # values near the top left
                if (row+col < cop_coeff_2*rows -1):
                    red_mass = red_mass + copula[row][col]
                    top_left = top_left + copula[row][col]
                    
                # values near the bottom right
                elif (row+col > cop_coeff_3*rows -1):
                    red_mass = red_mass + copula[row][col]
                    bottom_right = bottom_right + copula[row][col]
            
            # values in the other diagonal
            else:
                # values near the top right and bottom left corner
                if (row+col >= cop_coeff_2*rows -1) & (row+col <= cop_coeff_3*rows -1):                    
                    # values near the top right
                    if row <= rows / 2:
                        blue_mass = blue_mass + copula[row][col]
                        top_right = top_right + copula[row][col]
                    
                    # values near the bottom left
                    elif row > rows / 2:
                        blue_mass = blue_mass + copula[row][col]
                        bottom_left = bottom_left + copula[row][col]
    
    
    indicator = (red_mass+1e-9) / (blue_mass+1e-9)
    
    if indicator > indicator_cutoff:
        if (top_left / bottom_right) > indicator_cutoff:
            dependence = "positive_upper_tail"
        elif (bottom_right / top_left) > indicator_cutoff:
            dependence = "positive_lower_tail"
        else:
            dependence = "positive_upper_lower_tail"
        
    elif indicator < (1/indicator_cutoff):
        if (top_right / bottom_left) > indicator_cutoff:
            dependence = "negative_upper_tail"
        elif (bottom_left / top_right) > indicator_cutoff:
            dependence = "negative_lower_tail"
        else:
            dependence = "negative_upper_lower_tail"
            
    else:
        dependence = "no_dependence"
            
    
    return dependence, indicator


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


def find_reactants_products(cobra_model, reactions_ids=[]):
    
    #reactions_ids =  [ reaction.id for reaction in cobra_model.reactions ]
    
    reactants_list_all_reactions = []
    products_list_all_reactions = []
    reversibility_list_all_reactions = []
    
    for reaction in reactions_ids:
        reaction_id = reaction.id
        if reaction_id.endswith("_rev"):
            reaction_id = reaction[:-4]
        
        reactants_list_single_reaction = []
        products_list_single_reaction = []

        reaction_information = cobra_model.reactions.get_by_id(reaction_id)
        
        reactants = reaction_information.reactants
        products = reaction_information.products
        
        reversibility = reaction_information.reversibility
        reversibility_list_all_reactions.append(reversibility)
        
        for reactant in reactants:
            reactant = str(reactant)
            reactants_list_single_reaction.append(reactant)
                
        for product in products:
            product = str(product)
            products_list_single_reaction.append(product)         
                
        reactants_list_all_reactions.append(reactants_list_single_reaction)
        products_list_all_reactions.append(products_list_single_reaction)
        
    return reversibility_list_all_reactions, reactants_list_all_reactions, products_list_all_reactions


def sharing_metabolites(reactions_ids=[], reversibility_list_all_reactions=[], reactants_list_all_reactions=[], products_list_all_reactions=[], reaction_a="", reaction_b=""):
        
    #reactions_ids =  [ reaction.id for reaction in cobra_model.reactions ]

    reaction_a_index = reactions_ids.index(reaction_a)
    reaction_b_index = reactions_ids.index(reaction_b)
    
    reactants_a = reactants_list_all_reactions[reaction_a_index]
    reactants_b = reactants_list_all_reactions[reaction_b_index]
    
    products_a = products_list_all_reactions[reaction_a_index]
    products_b = products_list_all_reactions[reaction_b_index]
    
    if reaction_a_index == reaction_b_index:
            return True
        
    if (reversibility_list_all_reactions[reaction_a_index] == True) or (reversibility_list_all_reactions[reaction_b_index] == True):

        for reactant_a in reactants_a:
            if (reactant_a in reactants_b) or (reactant_a in products_b):
                return True
            
        for product_a in products_a:
            if (product_a in reactants_b) or (product_a in products_b):
                return True
                    
    # not reversible reactions
    elif (reversibility_list_all_reactions[reaction_a_index] == False) and (reversibility_list_all_reactions[reaction_b_index] == False):
        
        for reactant_a in reactants_a:
            if (reactant_a in products_b):
                return True
                
        # keep searching for sharing metabolites from products of A
        for product_a in products_a:
            if product_a in reactants_b:
                return True
            

def sharing_metabolites_square_matrix(reactions_ids=[], reversibility_list_all_reactions=[], reactants_list_all_reactions=[], products_list_all_reactions=[]):
    reactions_count = len(reactions_ids)
    boolean_sharing_metabolites_matrix = np.full((reactions_count, reactions_count), None, dtype=bool)

    for reaction_a in reactions_ids:
        for reaction_b in reactions_ids:
            reaction_a_id = reaction_a.id
            reaction_b_id = reaction_b.id
            
            reaction_a_id_index = reactions_ids.index(reaction_a_id)
            reaction_b_id_index = reactions_ids.index(reaction_b_id)
        
            sharing = sharing_metabolites(reactions_ids, reversibility_list_all_reactions, 
                                reactants_list_all_reactions, products_list_all_reactions,
                                reaction_a_id, reaction_b_id)
            
            boolean_sharing_metabolites_matrix[reaction_a_id_index, reaction_b_id_index] = sharing
            boolean_sharing_metabolites_matrix[reaction_b_id_index, reaction_a_id_index] = sharing

    return boolean_sharing_metabolites_matrix
