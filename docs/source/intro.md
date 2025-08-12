
# Introduction

## Introduction to dingo

[dingo](https://github.com/GeomScale/dingo) is a Python package that analyzes metabolic networks. It relies on high dimensional sampling with Markov Chain Monte Carlo (MCMC) methods and fast optimization methods to analyze the possible states of a metabolic network. To perform MCMC sampling, dingo relies on the C++ library volesti, which provides several algorithms for sampling convex polytopes. dingo also performs two standard methods to analyze the flux space of a metabolic network, namely Flux Balance Analysis and Flux Variability Analysis.


## Introduction to dingo-stats

[dingo-stats](https://github.com/SotirisTouliopoulos/dingo-stats) is a Python package for the statistical analysis of flux sampling results. It works given any flux sampling dataset (in a desired format) and thus is not restricted only to dingo samples. Functions from `dingo-stats` are used for identifying and removing loopy reactions (as a quality control step) and estimating their effect on thermodynamic feasibility of the samples. Additionally users can extract pathway information from the KEGG database and compare different sampling datasets with differential flux and pathways enrichment analyses. Moreover, users can calcate pairwise linear correlations and non-linear dependencies using copulas between reaction fluxes, perform and evaluate clustering results for metabolic pathways/operons prediction and proceed to graph analysis to reveal patterns from node (reaction) centrality metrics with flux distribution stats or reaction essentialities.