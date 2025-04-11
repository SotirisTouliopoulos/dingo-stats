
from cobra.flux_analysis import flux_variability_analysis
from cobra.io import read_sbml_model


ec_model = read_sbml_model('../../data/Reconstructions_data/e_coli_core.xml')
reaction_list = ec_model.reactions


fva = flux_variability_analysis(ec_model, loopless=False, fraction_of_optimum=0)
#print(fva)
fva_dict = fva.apply(lambda row: (row['minimum'], row['maximum']), axis=1).to_dict()
#print(fva_dict)

fva_loopless = flux_variability_analysis(ec_model, loopless=True, fraction_of_optimum=0)
print(fva_loopless)
fva_loopless_dict = fva.apply(lambda row: (row['minimum'], row['maximum']), axis=1).to_dict()


print(ec_model.reactions.get_by_id(reaction).bounds for reaction in ec_model.reactions)

for reaction in reaction_list:
    model_reaction = ec_model.reactions.get_by_id(reaction.id)
    model_reaction.bounds = (fva_loopless_dict[reaction.id][0], fva_loopless_dict[reaction.id][1])

print(ec_model.reactions.get_by_id(reaction).bounds for reaction in ec_model.reactions)



