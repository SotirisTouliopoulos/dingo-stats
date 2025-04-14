
from cobra.io import load_matlab_model
from dingo import MetabolicNetwork


roseburia_cobra_model = load_matlab_model("/home/touliopoulos/project/erasmus/erasmus_2025_project/data/Reconstructions_data/Roseburia_intestinalis_L1_82_fixed.mat")
roseburia_dingo_model = MetabolicNetwork.from_mat("/home/touliopoulos/project/erasmus/erasmus_2025_project/data/Reconstructions_data/Roseburia_intestinalis_L1_82_fixed.mat")


metabolites = roseburia_dingo_model.metabolites
# "ac[e]" 
# "lac_l[e]" 
# "lac_d[e]" 
# "glc_d[e]" 
# "fru[e]" 
# "gal[e]" 
# "sucr[e]" 
# "malt[e]"


metabolites_matches = [m for m in metabolites if "malt[e]" == m.lower()]


def find_metabolite_in_pathway_groups(cobra_model, desired_metabolite):
    for group in cobra_model.groups:
        for reaction in group.members:
            for metabolite in reaction.metabolites:
                if str(metabolite).lower() == desired_metabolite:
                    print(group, reaction)


find_metabolite_in_pathway_groups(roseburia_cobra_model, "ac[e]")
find_metabolite_in_pathway_groups(roseburia_cobra_model, "lac_d[e]")
find_metabolite_in_pathway_groups(roseburia_cobra_model, "lac_l[e]")

find_metabolite_in_pathway_groups(roseburia_cobra_model, "glc_d[e]")
find_metabolite_in_pathway_groups(roseburia_cobra_model, "fru[e]")
find_metabolite_in_pathway_groups(roseburia_cobra_model, "gal[e]")
find_metabolite_in_pathway_groups(roseburia_cobra_model, "sucr[e]")
find_metabolite_in_pathway_groups(roseburia_cobra_model, "malt[e]")

