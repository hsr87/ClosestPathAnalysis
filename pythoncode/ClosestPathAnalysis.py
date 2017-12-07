# networkx package should be installed
import networkx as nx

def closestPathAnalysis(set_ds_input, lst_ta_input, grp_input):
    for diseasename in set_ds_input:
        f_ds_score = open(
            "../data/result_file/tr_ds_dist/Scores_" + diseasename + ".txt",
            'w')
        str_header = "Target" + "\t" + "Distance\n"
        f_ds_score.write(str_header)
        for drugtarget in lst_ta_input:
            try:
                int_target_dist = nx.shortest_path_length(grp_input, source=drugtarget, target=diseasename)
            except nx.exception.NetworkXNoPath:
                int_target_dist = 10000
                pass
            str_line = drugtarget + "\t" + str(int_target_dist) + "\n"
            f_ds_score.write(str_line)
        f_ds_score.close()
        print (diseasename + " is done!")

if __name__ == '__main__':
    # Read PPI
    f_ppi = open("../data/input_file/BioGRID_PPI.tsv", 'r')
    lst_relations = f_ppi.readlines()
    f_ppi.close()
    print('Read PPI sucessfully')
    lst_relations.pop(0)
    SPLIT = list(map(lambda relation: relation.split("\t"), lst_relations))
    lst_left = [ele[0] for ele in SPLIT]
    lst_right = [ele[1] for ele in SPLIT]
    lst_right = list(map(lambda s: s.rstrip(), lst_right))
    lst_rel = ['ppi'] * len(lst_left)
    set_left = set(lst_left)
    set_right = set(lst_right)
    set_nodes_ppi = set_left.union(set_right)

    # Read a disease list
    f_ds = open("../data/input_file/CTD_genes_diseases_direct.tsv", 'r')
    lst_diseases = f_ds.readlines()
    f_ds.close()
    print('Read DS sucessfully')
    lst_diseases.pop(0)
    SPLIT = list(map(lambda disease: disease.split("\t"), lst_diseases))
    lst_dg = [ele[0] for ele in SPLIT]
    lst_ds = [ele[2] for ele in SPLIT]
    lst_ds = list(map(lambda s: s.rstrip(), lst_ds))
    set_ds = set(lst_ds)
    set_nodes = set_nodes_ppi.union(set_ds)

    # Read a drug list
    f_dr = open("../data/input_file/target_lists.tsv", 'r')
    lst_drugs = f_dr.readlines()
    f_dr.close()
    print('Read CP sucessfully')
    SPLIT = list(map(lambda drug: drug.split("\t"), lst_drugs))
    lst_ta_notfiltered = [ele[0] for ele in SPLIT]
    lst_ta_notfiltered = list(map(lambda s: s.rstrip(), lst_ta_notfiltered))

    lst_ta = []
    for nIndex in range(len(lst_ta_notfiltered)):
        if (lst_ta_notfiltered[nIndex] in set_nodes):
            lst_ta.append(lst_ta_notfiltered[nIndex])

    # Network construction
    grp_ppi = nx.DiGraph()
    grp_ppi.add_nodes_from(set_nodes)
    lst_edges = []
    for nIndex in range(len(lst_rel)):
        lst_edges.append((lst_left[nIndex], lst_right[nIndex]))
        lst_edges.append((lst_right[nIndex], lst_left[nIndex]))
    for nIndex in range(len(lst_dg)):
        lst_edges.append((lst_dg[nIndex], lst_ds[nIndex]))
    grp_ppi.add_edges_from(lst_edges)

    del lst_edges,  lst_ta_notfiltered, SPLIT, lst_drugs, lst_diseases, lst_relations, lst_rel

    closestPathAnalysis(set_ds, lst_ta, grp_ppi)


