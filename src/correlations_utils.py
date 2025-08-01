
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import math
from scipy.spatial.distance import jensenshannon
from scipy.stats import spearmanr
from typing import Dict, Tuple, List
from cobra import Model, Reaction
from numpy.typing import NDArray


BIGG_COFACTORS = ['atp_c0', 'atp_c', 'adp_c', 'adp_c0',
                  'atp_c0', 'atp_c', 'adp_c', 'adp_c0',
                  'udp_c0', 'udp_c', 'ump_c0', 'ump_c',
                  'amp_c', 'amp_c0',
                  'gdp_c0', 'gdp_c', 'gtp_c0', 'gtp_c',
                  'accoa_c', 'accoa_c0', 'coa_c', 'coa_c0',  # acetyl-CoA
                  'q8_c0', 'q8_c', 'q8h2_c', 'q8h2_c0', 'mqn8_c', 'mqn8_c0', 'mql8_c', 'mql8_c0', 'q8h2_c', 'q8h2_c0',
                  'actp_c0', 'actp_c',
                  'h2o_c', 'h2o_c0', 'h2o_e', 'h2o[e]',
                  'pi_e', 'pi[e]', 'pi_c', 'pi_c0', 'ppi_c0', 'ppi_c',
                  'pep_c', 'pep_c0',
                  'h_c', 'h_c0', 'h_e', 'h[e]',
                  'o2_c', 'o2_c0', 'o2_e', 'o2[e]',
                  'co2_c', 'co2_c0', 'co2_e', 'co2[e]',
                  'nadp_c', 'nadp_c0', 'nadph_c', 'nadph_c0', 'nad_c', 'nad_c0', 'nadh_c', 'nadh_c0',
                  'nadp_e', 'nadp[e]', 'nadph_e', 'nadph_c0', 'nad_e', 'nad[e]', 'nadh_e', 'nadh[e]',
                  'fadh2_c', 'fadh2_c0', 'fad_c', 'fad_c0',
                  'nh4_c', 'nh4_c0', 'nh4_e', 'nh4[e]',
                  'pyr_c0', 'pyr_c'
                  ]


BIGG_BUILDING_BLOCKS = ['ala_L_c0', 'asp_L_c0', ' gln_L_c0', 'glu_L_c0',
                        'glu_L_c0', 'ser_L_c0', 'trp_L_c0', 'met_L_c0',
                        'lys_L_c0', 'cyst_L_c0']


MODELSEED_COFACTORS = [
    "cpd00001_c0",  # h2o
    "cpd00002_c0",  # atp
    "cpd00003_c0",  # nad
    "cpd00004_c0",
    "cpd00005_c0",
    "cpd00006_c0",  # nadp
    "cpd00007_c0",
    "cpd00008_c0",  # adp
    "cpd00009_c0",  # HZ added
    "cpd00010_c0",  # CoA
    "cpd00011_c0",  # co2
    "cpd00012_c0",  # ppi
    "cpd00013_c0",  # NH3
    "cpd00014_c0",
    "cpd00015_c0",  # fad
    "cpd00018_c0",  # amp-like
    "cpd00020_c0",  # pyruvate
    "cpd00022_c0",
    "cpd00031_c0",  # gdp-like
    "cpd00038_c0",  # gtp
    "cpd00056_c0",  # ttp
    "cpd00061_c0",  # pep
    "cpd00067_c0",  # H+
    "cpd15353_c0",
    "cpd15499_c0",
    "cpd15561_c0",
    "cpd00097_c0",
    "cpd00982_c0",
    "cpd01270_c0",
    "cpd00052_c0",
    "cpd00062_c0",
    "cpd00068_c0",
    "cpd00115_c0",
    "cpd00241_c0",
    "cpd00356_c0",
    "cpd00357_c0",
    "cpd00358_c0",
    "cpd00530_c0",
    "cpd00977_c0",
    "cpd01775_c0"
]


def correlated_reactions(
    steady_states: NDArray[np.float64], 
    reactions: List = [], 
    linear_coeff: str = "pearson",
    include_non_linear: bool = False, 
    boolean_sharing_metabolites_matrix: NDArray[np.float64] = None,
    linear_corr_cutoff: float = 0.30, 
    indicator_cutoff: float = 1.2, 
    jensenshannon_cutoff: float = 0.1, 
    std_cutoff: float = 1e-3,
    cells: int = 4,
    cop_coeff: float = 0.2,
    lower_triangle: bool = True, 
    verbose: bool = False
) -> Tuple:
    """
    Function that calculate pairwise linear correlations and non-linear copula dependencies given a flux sampling array
    User can choose the preffered coefficient to calculate pairwise linear correlations (pearson or spearman). 
    Moreover, he can filter correlations not meeting a cutoff value and choose whether to proceed with calculation of
    non-linear dependencies using copulas.

    Keyword arguments:
    steady_states (NDArray[np.float64]) -- A numpy array of the generated steady states fluxes
    reactions (List) -- A list with the reactions IDs (must be in accordance with the rows of the steady states)
    linear_coeff (float) -- A string declaring the linear coefficient to be selected. Available options: "pearson", "spearman"
    include_non_linear (bool) -- A boolean variable that if True, takes into account and calculates non-linear correlations
    boolean_sharing_metabolites_matrix (NDArray[np.float64]) -- A boolean symmetric numpy 2D array with True/False based on the presense of shared metabolites between reactions
    linear_corr_cutoff (float) -- A cutoff to filter (remove) linear correlations (not greater than the cutoff) based on the pearson or spearmanr coefficient
    indicator_cutoff (float) -- A cutoff to classify non-linear correlations as positive, negative or non-significant
    jensenshannon_cutoff (float) -- A cutoff to filter (remove) non-linear correlations (not greater than the cutoff) based on the Jensen-Shannon metric
    std_cutoff (float) -- A cutoff to avoid computing the copula between 2 fluxes with almost fixed values
    cells (int) -- Number of cells to compute the copula
    cop_coeff (float) -- A value that narrows or widens the width of the copula's diagonal (use lower values to capture extreme tail dependences)
    lower_triangle (bool) -- A boolean variable that if True returns only the lower triangular matrix
    verbose (bool) -- A boolean variable that if True additional information is printed as an output.

    Returns:
    if include_non_linear is set to False:
    Tuple[NDArray[np.float64], Dict]
        linear_correlation_matrix -- A correlations matrix filtered based on the given cutoffs that includes only linear correlations
        correlations_dictionary -- A dictionary containing unique reaction pairs and their corresponding correlation values

    if include_non_linear is set to True:
    Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], Dict]
        linear_correlation_matrix -- A correlations matrix filtered based on the given cutoffs that includes only linear correlations
        correlations_dictionary -- A dictionary containing unique reaction pairs and their corresponding correlation values
        mixed_correlation_matrix -- A correlations matrix filtered based on the given cutoffs that includes both linear and non-linear correlations
        non_linear_correlation_matrix -- A correlations matrix filtered based on the given cutoffs that includes only non-linear correlations
    """

    if indicator_cutoff < 1:
        raise Exception("Indicator cutoff must be at least equal to 1")

    if cop_coeff > 0.4 or cop_coeff < 0.2:
        raise Exception("Input value to cop_coeff parameter must be between 0.2 and 0.4")

    # calculate coefficients to access red and blue copula mass
    cop_coeff_1 = cop_coeff
    cop_coeff_2 = 1 - cop_coeff
    cop_coeff_3 = 1 + cop_coeff

    # compute correlation matrix
    if linear_coeff == "pearson":
        corr_matrix = np.corrcoef(steady_states, rowvar=True)
    elif linear_coeff == "spearman":
        corr_matrix, _ = spearmanr(steady_states, axis=1)
    else:
        raise Exception("Input value to linear_coeff parameter is not valid. Choose between pearson or spearman")
    
    # fill diagonal with 1
    np.fill_diagonal(corr_matrix, 1)
    # replace not assigned values with 0
    corr_matrix[np.isnan(corr_matrix)] = 0
    # keep only the lower triangle (with the diagonal) to reduce computational time
    corr_matrix[np.triu_indices(corr_matrix.shape[0], 1)] = np.nan

    # create a copy of correlation matrix to replace/filter values
    linear_correlation_matrix = corr_matrix.copy()

    # find indices of correlation matrix where correlation does not occur
    no_corr_indices = np.argwhere((linear_correlation_matrix < linear_corr_cutoff) & (linear_correlation_matrix > -linear_corr_cutoff))
    # find indices of correlation matrix where correlation does occur
    corr_indices = np.argwhere((linear_correlation_matrix > linear_corr_cutoff) | (linear_correlation_matrix < -linear_corr_cutoff))

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

            # define a hypothetical uniform copula to compare with the previously computed copula
            reference_copula = np.full((cells, cells), (1 / (cells * cells)))
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
                print("Completed the process of", i+1 ,"from",no_corr_indices.shape[0], "copulas")


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


def plot_correlation_matrix(
    correlation_matrix: NDArray[np.float64], 
    reactions: List = [], 
    label_font_size: int = 5,
    width: int = 900,
    height: int = 900
):
    """
    Function that plots a correlation matrix created with the `correlated_reactions` function.

    Keyword arguments:
    correlation_matrix (2D array) -- The correlation matrix to plot.
    reactions (list) -- Optional list of reaction labels for axes.
    label_font_size (int) -- Font size for axis labels (default: 5).
    width (int) -- Width of the plot
    height (int) -- Height of the plot
    """

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
        zmin=-1,
        zmax=1,
        x=reactions,
        y=reactions,
        origin="upper"
    )

    fig.update_layout(
        xaxis=dict(tickfont=dict(size=label_font_size)),
        yaxis=dict(tickfont=dict(size=label_font_size)),
        width=width,
        height=height,
        plot_bgcolor="rgba(0,0,0,0)"
    )

    fig.update_traces(xgap=1, ygap=1, hoverongaps=False)
    fig.show()


def compute_copula(
    flux1: NDArray[np.float64], 
    flux2: NDArray[np.float64],
    n: int
) -> NDArray[np.float64]:
    """
    Function that estimates the copula between two fluxes

    Keyword arguments:
    flux1 (NDArray[np.float64]) -- A vector that contains the measurements of the first reaxtion flux
    flux2 (NDArray[np.float64]) -- A vector that contains the measurements of the second reaxtion flux
    n (int) -- The number of cells
    
    Returns:
    copula (NDArray[np.float64]) -- The computed copula matrix 
    """

    N = flux1.size
    copula = np.zeros([n,n], dtype=float)

    I1 = np.argsort(flux1)
    I2 = np.argsort(flux2)

    grouped_flux1 = np.zeros(N)
    grouped_flux2 = np.zeros(N)

    for j in range(n):
        rng = range((j * math.floor(N / n)),((j + 1) * math.floor(N / n)))
        grouped_flux1[I1[rng]] = j
        grouped_flux2[I2[rng]] = j

    for i in range(n):
        for j in range(n):
            copula[i,j] = sum((grouped_flux1 == i) * ( grouped_flux2 == j))

    copula = copula / N
    return copula


def copula_tail_dependence(
    copula: NDArray[np.float64], 
    cop_coeff_1: float, 
    cop_coeff_2: float, 
    cop_coeff_3: float, 
    indicator_cutoff: float
) -> Tuple[str, float]:
    """
    Function that given a copula and parameters that define the diagonals width, classifies the tail-dependence between a reaction pair.

    Suppose we have a 5*5 copula (for better vizualization, this works with any dimensions) looking like:
    
    [
    [0.001      0.03766667 0.05       0.05433333 0.057     ]
    [0.015      0.04733333 0.04333333 0.048      0.04633333]
    [0.03       0.05       0.046      0.036      0.038     ]
    [0.05766667 0.03833333 0.035      0.034      0.035     ]
    [0.09633333 0.02666667 0.02566667 0.02766667 0.02366667] 
    ]

    - This copula array for our algorithm has 4 separate edges/areas:
        - Top left, for example cells close to: [0,0], [0,1], [1,0] matrix position
        - Bottom right, for example cells close to: [4,4], [4,3], [3,4] matrix position
        - Top right, for example cells close to: [4,0], [3,0], [4,1] matrix position
        - Bottom left, for example cells close to: [0,4], [0,3], [1,4] matrix position

    - This copula array for our algorithm has 2 diagonals: 
        - 1st diagonal: includes top left and bottom right areas
        - 2nd diaognal: includes top right and bottom left areas
        
    This function parses the 2 copula's diagonals and picks values from the 4 copula's edges. 
    One diagonal (the 1st) corresponds to the mass probability of 2 reactions working together at the same time 
    (when one is active the other is active too) and the other diagonal (the 2nd) corresponds to the mass propability 
    of 2 reactions not working together at the same time (when one is active the other is inactive).
    
    Then, by dividing the sum of propability mass of the 1st to the 2nd diagonal the indicator value is calculated. 
    This value informs us on which diagonal has a higher mass concentration. If the indicator value is above a significance cutoff, 
    this gives a sign to the dependence (positive for indicator value above given threshold or negative for indicator value 
    below the reciprocal of the given threshold).

    Division of the corresponding edges (e.g. for the 1st diagonal: top left with the bottom right edge) shows whether 
    higher mass concentration appears in one or both edges. This helps with classification of the copula dependence.

    We have different classification types based on:
    - sign: `positive` or `negative`, as explained above
    - tail(s): show which tail or tails have significant higher mass concentration based on the sign:
        - `positive_upper_tail` for higher mass in the top left area
        - `positive_lower_tail` for higher mass in the bottom right area
        - `positive_upper_lower_tail` for higher mass in both bottom right and top left areas
        
        - `negative_upper_tail` for higher mass in the top right area
        - `negative_lower_tail` for higher mass in the bottom left area
        - `negative_upper_lower_tail` for higher mass in both top right and bottom left areas

    Keyword arguments:
    cop_coeff (float) -- Parameters that define the width of the copula's diagonal (leading to more wide or narrow diagonals).
    indicator_cutoff (float) -- Cutoff to filter and reveal the sign (positive/negative) between a dependence. Must be greater than 1
    
    Returns:
    Tuple[str, float]
        dependence (str) -- The classification for the tail-dependence
        indicator (float) -- The value of the indicator
    """
    
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
            if ((row - col >= -cop_coeff_1 * rows) & (row - col <= cop_coeff_1 * rows)): 
                # values near the top left
                if (row + col < cop_coeff_2 * rows - 1):
                    red_mass = red_mass + copula[row][col]
                    top_left = top_left + copula[row][col]

                # values near the bottom right
                elif (row + col > cop_coeff_3 * rows -1):
                    red_mass = red_mass + copula[row][col]
                    bottom_right = bottom_right + copula[row][col]

            # values in the other diagonal
            else:
                # values near the top right and bottom left corner
                if (row + col >= cop_coeff_2 * rows - 1) & (row + col <= cop_coeff_3 * rows - 1):
                    # values near the top right
                    if row < rows / 2:
                        blue_mass = blue_mass + copula[row][col]
                        top_right = top_right + copula[row][col]

                    # values near the bottom left
                    elif row >= rows / 2:
                        blue_mass = blue_mass + copula[row][col]
                        bottom_left = bottom_left + copula[row][col]


    indicator = (red_mass + 1e-9) / (blue_mass + 1e-9)

    if indicator > indicator_cutoff:
        if (top_left / bottom_right) > indicator_cutoff:
            dependence = "positive_upper_tail"
        elif (bottom_right / top_left) > indicator_cutoff:
            dependence = "positive_lower_tail"
        else:
            dependence = "positive_upper_lower_tail"

    elif indicator < (1 / indicator_cutoff):
        if (top_right / bottom_left) > indicator_cutoff:
            dependence = "negative_upper_tail"
        elif (bottom_left / top_right) > indicator_cutoff:
            dependence = "negative_lower_tail"
        else:
            dependence = "negative_upper_lower_tail"

    else:
        dependence = "no_dependence"


    return dependence, indicator


def plot_copula(
    data_flux1: NDArray[np.float64], 
    data_flux2: NDArray[np.float64], 
    n: int = 5, 
    width: int = 900 , 
    height: int = 600, 
    export_format: str = "svg"
):
    """
    Function that plots the copula between two fluxes

    Keyword arguments:
    data_flux1 -- A list that contains: (i) the vector of the measurements of the first reaction,
                                      (ii) the name of the first reaction
    data_flux2 -- A list that contains: (i) the vector of the measurements of the second reaction,
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
        title = 'Copula between ' + data_flux1[1] + ' and ' + data_flux2[1],
        scene = dict(
            xaxis_title=data_flux1[1],
            yaxis_title=data_flux2[1],
            zaxis_title="prob, mass"
        ),
        margin=dict(r=30, b=30, l=30, t=50)
    )

    fig.layout.template = None
    fig.show()

    camera = dict(
        up     = dict(x=0, y=0, z=1),
        center = dict(x=0, y=0, z=0),
        eye    = dict(x=1.25, y=1.25, z=1.25)
    )

    fig.update_layout(scene_camera=camera)
    fig.to_image(format = export_format, engine="kaleido")


def split_forward_reverse(
    steady_states: NDArray[np.float64], 
    reactions: List = []
) -> Tuple[NDArray[np.float64], List]:
    """
    Function that given a flux sampling (steady states) array with reactions as rows, splits 
    all bidirectional reactions having at least 1 positive and 1 negative flux value into separate 
    forward and reverse reactions.
    
    Keyword arguments:
    steady_states (NDArray[np.float64]) -- Flux Sampling array with only forward reactions
    reactions (List) -- List with reactions of interest to split and create separate forward-reverse reactions 
                        (pathway of interest, not all reactions will be split)
                        
    Return: 
    extended_steady_states (NDArray[np.float64]) -- Extended Flux Sampling array with separate forward and reverse fluxes in bidirectional reactions
    extended_reactions (List) -- Extended List including the names of the added reverse reactions 
    """

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


def find_reactants_products(
    cobra_model: Model, 
    reactions_ids: List = []
) -> Tuple[List, List, List]:
    """
    Function that identifies and saves in separate lists the reactants, products and directionality status 
    of each reaction in a given metabolic model. Cofactors do not count as reactants/products.
    
    Keyword arguments:
    cobra_model (Model) -- Metabolic Model in cobra format
    reactions_ids (List) -- List of the reactions of interest
    
    Returns:
    Tuple[List, List, List]
        reversibility_list_all_reactions (List) -- List containing the reversibility status of the given reaction IDs
        reactants_list_all_reactions (List) -- List containing the reactants of the given reaction IDs
        products_list_all_reactions (List) -- List containing the products of the given reaction IDs
    """

    reactants_list_all_reactions = []
    products_list_all_reactions = []
    reversibility_list_all_reactions = []

    for reaction in reactions_ids:
        if reaction.endswith("_rev"):
            reaction_id = reaction[:-4]
        else:
            reaction_id = reaction

        reactants_list_single_reaction = []
        products_list_single_reaction = []

        reaction_information = cobra_model.reactions.get_by_id(reaction_id)

        reactants = reaction_information.reactants
        products = reaction_information.products

        reversibility = reaction_information.reversibility
        reversibility_list_all_reactions.append(reversibility)

        for reactant in reactants:
            reactant = str(reactant)
            if (reactant not in BIGG_COFACTORS) and (reactant not in MODELSEED_COFACTORS):
                reactants_list_single_reaction.append(reactant)

        for product in products:
            product = str(product)
            if (product not in BIGG_COFACTORS) and (product not in MODELSEED_COFACTORS):
                products_list_single_reaction.append(product)

        reactants_list_all_reactions.append(reactants_list_single_reaction)
        products_list_all_reactions.append(products_list_single_reaction)

    return reversibility_list_all_reactions, reactants_list_all_reactions, products_list_all_reactions


def sharing_metabolites(
    reactions_ids: List = [], 
    reversibility_list_all_reactions: List = [], 
    reactants_list_all_reactions: List = [], 
    products_list_all_reactions: List = [], 
    reaction_a: str = "", 
    reaction_b: str = ""
) -> bool:
    """
    Function that compares the reactants and products of 2 reactions and returns True 
    when a common metabolite (reactant/product and not cofactor) exists between them
    
    Keyword arguments:
    reactions_ids (List) -- List containing IDs from a set of reactions of interest
    reversibility_list_all_reactions (List) -- List containing the reversibility status from a set of reactions 
                                               (Reactions must be the same and in accordance with their order in `reactions_ids` list)
    reactants_list_all_reactions (List) -- List containing the reactants from a set of reactions
                                           (Reactions must be the same and in accordance with their order in `reactions_ids` list)
    products_list_all_reactions (List) -- List containing the products from a set of reactions
                                          (Reactions must be the same and in accordance with their order in `reactions_ids` list)
    reaction_a (str) -- ID of reaction 1
    reaction_b (str) -- ID of reaction 2
    
    Returns:
    bool -- True/False based on the presence of a common metabolite between the 2 reactions
    """

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


def sharing_metabolites_square_matrix(
    reactions_ids: List = [],
    reversibility_list_all_reactions: List = [],
    reactants_list_all_reactions: List = [],
    products_list_all_reactions: List = []
) -> NDArray[np.float64]:    
    """
    Function that given the lists with reactants, products and reversibility information,
    creates a square boolean matrix with True values representing reactions sharing 
    a common metabolite, as implemented in `sharing_metabolites` function.
    
    Keyword arguments:
    reactions_ids (List) -- List containing IDs from a set of reactions of interest
    reversibility_list_all_reactions (List) -- List containing the reversibility status from a set of reactions 
                                               (Reactions must be the same and in accordance with their order in `reactions_ids` list)
    
    reactants_list_all_reactions (List) -- List containing the reactants from a set of reactions
                                           (Reactions must be the same and in accordance with their order in `reactions_ids` list)
    
    products_list_all_reactions (List) -- List containing the products from a set of reactions
                                          (Reactions must be the same and in accordance with their order in `reactions_ids` list)
    
    Returns:
    boolean_sharing_metabolites_matrix (NDArray[np.float64]) -- Boolean Square Matrix with information on reactions sharing metabolites
    """

    reactions_count = len(reactions_ids)
    boolean_sharing_metabolites_matrix = np.full((reactions_count, reactions_count), None, dtype=bool)

    for reaction_a in reactions_ids:
        for reaction_b in reactions_ids:
            #reaction_a_id = reaction_a.id
            #reaction_b_id = reaction_b.id

            reaction_a_id_index = reactions_ids.index(reaction_a)
            reaction_b_id_index = reactions_ids.index(reaction_b)
            
            sharing = sharing_metabolites(
                reactions_ids, reversibility_list_all_reactions,
                reactants_list_all_reactions, products_list_all_reactions,
                reaction_a, reaction_b
            )

            boolean_sharing_metabolites_matrix[reaction_a_id_index, reaction_b_id_index] = sharing
            boolean_sharing_metabolites_matrix[reaction_b_id_index, reaction_a_id_index] = sharing

    return boolean_sharing_metabolites_matrix