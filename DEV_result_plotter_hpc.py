
import matplotlib.pyplot as plt

def sort_results(index_list, list_): #puts in energy order
    return [list_[i] for i in index_list]

def main(results_list, main_directory):
    """Produces results text file and graphs for the collection of simulations. """

    #extracting data from results dict from atom_track.py
    energies = [float(results_dict['energy']) for results_dict in results_list]
    sorted_energies = sorted(energies, key = float)
    index_list = [energies.index(energy) for energy in sorted_energies]
    
    energies = sorted_energies
    avg_pens = sort_results(index_list,[float(results_dict['avg_pen'][0]) for results_dict in results_list])
    pen_errs = sort_results(index_list,[float(results_dict['avg_pen'][1]) for results_dict in results_list])
    diamond_avg_pens = sort_results(index_list,[float(results_dict['diamond_avg_pen'][0]) for results_dict in results_list])
    diamond_pen_errs = sort_results(index_list,[float(results_dict['diamond_avg_pen'][1]) for results_dict in results_list])
    max_pens = sort_results(index_list,[results_dict['max_pen'] for results_dict in results_list])
    diamond_max_pens = sort_results(index_list,[results_dict['diamond_max_pen'] for results_dict in results_list])
    repeats = sort_results(index_list,[results_dict["repeats"] for results_dict in results_list])
    bulk_diamond_percs = sort_results(index_list,[float(results_dict['diamond_atoms_lost'][0]) for results_dict in results_list])
    bulk_diamond_percs_errs = sort_results(index_list,[float(results_dict['diamond_atoms_lost'][1]) for results_dict in results_list])
    bulk_graphene_percs = sort_results(index_list,[float(results_dict['graphene_atoms_lost'][0]) for results_dict in results_list])
    bulk_graphene_percs_errs = sort_results(index_list,[float(results_dict['graphene_atoms_lost'][1]) for results_dict in results_list])

    no_of_sheets = float(results_list[0]['no_of_sheets'])
    atom_type = results_list[0]["atom_type"]
    diamond_type = results_list[0]["diamond_type"]

    pens_list = sort_results(index_list,[results_dict['pens'] for results_dict in results_list])
    locations_list = sort_results(index_list,[results_dict['locations'] for results_dict in results_list])

    diamond_pen_percs = [pens[0]*100/repeats[pens_list.index(pens)] for pens in pens_list]
    graphene_pen_percs = [pens[1]*100/repeats[pens_list.index(pens)] for pens in pens_list]

    diamond_finish_percs = [finish[0]*100/repeats[locations_list.index(finish)] for finish in locations_list]
    graphene_finish_percs = [finish[1]*100/repeats[locations_list.index(finish)] for finish in locations_list]
    outside_finish_percs = [finish[2]*100/repeats[locations_list.index(finish)] for finish in locations_list]



 ################### Penetration Values Graph

    fig = plt.figure()

    plt.errorbar(x = energies, y = avg_pens, yerr=pen_errs, capsize= 5)
    plt.errorbar(x = energies, y = diamond_avg_pens, yerr=diamond_pen_errs, capsize= 5)
    plt.plot(energies, max_pens)

    plt.legend(["Maximum" , "Average", "Diamond Average"])
    
    plt.title("%s Graphene Sheets, %s Penetration Values"%(int(no_of_sheets), atom_type))
    plt.xlabel("Energy / eV")
    plt.ylabel("Penetration / Å")

    plt.savefig('%s/%s_%s_penetration_values.png'%(main_directory,atom_type, int(no_of_sheets)))
    plt.close(fig)
    
################### Final Location Graph

    fig = plt.figure()
    plt.plot(energies, diamond_finish_percs)
    plt.plot(energies, graphene_finish_percs)
    plt.plot(energies, outside_finish_percs)

    plt.legend(["Diamond" , "Graphene", "Outside"])

    plt.title("%s Graphene Sheets, %s Final Location"%(int(no_of_sheets),atom_type))
    plt.xlabel("Energy / eV")
    plt.ylabel("Percent")
    plt.ylim(-5,105)

    plt.savefig('%s/%s_%s_final_location.png'%(main_directory,atom_type, int(no_of_sheets)))
    plt.close(fig)

#################### Penetration percentages Graph

    fig = plt.figure()
    plt.plot(energies, diamond_pen_percs)
    if no_of_sheets != 0:
        plt.plot(energies, graphene_pen_percs)

    plt.legend(["Diamond" , "Graphene"])

    plt.title("%s Graphene Sheets, %s Penetration Percentages"%(int(no_of_sheets), atom_type))
    plt.xlabel("Energy / eV")
    plt.ylabel("Percent")
    plt.ylim(-5,105)

    plt.savefig('%s/%s_%s_penetration_perc.png'%(main_directory,atom_type, int(no_of_sheets)))
    plt.close(fig)

##################### Displaced Atoms graph

    fig = plt.figure()
    plt.errorbar(x = energies, y = bulk_diamond_percs, yerr=bulk_diamond_percs_errs, capsize= 5)
    plt.errorbar(x = energies, y = bulk_graphene_percs, yerr=bulk_graphene_percs_errs, capsize= 5)

    plt.legend(["Diamond" , "Graphene"])

    plt.title("%s Graphene Sheets, %s Displaced Atoms Percentages"%(int(no_of_sheets), atom_type))
    plt.xlabel("Energy / eV")
    plt.ylabel("Percent")


    plt.savefig('%s/%s_%s_displaced_atoms_perc.png'%(main_directory,atom_type,int(no_of_sheets)))
    plt.close(fig)

###################### Results text file - gives all data used in graphs

    data_str = "Data - %s Graphene Sheets, %s bombardment, %s diamond, %s repeats \n\n"%(int(no_of_sheets), atom_type, 
                                                                                        diamond_type, repeats[0])
    data_str += "\nGraph 1 - Penetration Values"
    data_str += '\n|{:12s} |{:12s} |{:12s} |{:12s} |{:12s} |{:12s} |'.format('Energy', "Avg. Pen.", "Error","Diamond", "Error", "Max Pen.")
    data_str += '\n|{:12s} |{:12s} |{:12s} |{:12s} |{:12s} |{:12s} |'.format('', "", "","Avg. Pen. ", "", "")

    for i in range(len(energies)):
        data_str += '\n|{:12.3f} |{:12.3f} |{:12.3f} |{:12.3f} |{:12.3f} |{:12.3f} |'.format(energies[i], avg_pens[i], pen_errs[i], diamond_avg_pens[i], diamond_pen_errs[i], max_pens[i])
    
    data_str += "\n\nGraph 2 - Finishing Location"
    data_str += '\n|{:12s} |{:12s} |{:12s} |{:12s} |'.format('Energy', "Diamond", "Graphene", "Outside")

    for i in range(len(energies)):
        data_str += '\n|{:12.3f} |{:12.3f} |{:12.3f} |{:12.3f} |'.format(energies[i], diamond_finish_percs[i], graphene_finish_percs[i], outside_finish_percs[i])

    data_str += "\n\nGraph 3 - Regions Penetrated"
    data_str += '\n|{:12s} |{:12s} |{:12s} |'.format('Energy', "Diamond", "Graphene")

    for i in range(len(energies)):
        data_str += '\n|{:12.3f} |{:12.3f} |{:12.3f} |'.format(energies[i], diamond_pen_percs[i], graphene_pen_percs[i])

    data_str += "\n\nGraph 4 - Displaced Atoms"
    data_str += '\n|{:12s} |{:12s} |{:12s} |{:12s} |{:12s} |'.format('Energy', "Diamond", "Error", "Graphene", "Error")

    for i in range(len(energies)):
        data_str += '\n|{:12.3f} |{:12.3f} |{:12.3f} |{:12.3f} |{:12.3f} |'.format(energies[i], bulk_diamond_percs[i], bulk_diamond_percs_errs[i], bulk_graphene_percs[i], bulk_graphene_percs_errs[i])

    try:
        with open("%s/%s_%s_%s_data.txt"%(main_directory, atom_type[0].lower(), int(no_of_sheets), diamond_type), 'w') as fp: 
            fp.write(data_str) 
    except FileNotFoundError:
        with open("%s/data.txt"%(main_directory), 'w') as fp: 
            fp.write(data_str) 
        







