import networkx as nx
import random
import numpy as np
import sys
import os
##sys.path.append('src/')


def generate_connected_graph_from_identity_score(identity_threshold):

    print("\n ------- Generating connected graph based on RNA sequence similarity ------- \n")
    os.chdir('src/')
    
    f_organism_name = open("Organism_list/Current_Organism_list","r")
    organism_name = []

    while(True):
        line = f_organism_name.readline()
        
        if line == "":
            break
        
        line = line.strip("\n")
        organism = line.replace(" ", "_")
        organism_name.append(organism)

    f_organism_name.close()


    for i in range(0, len(organism_name)):

        f_open_nodelist = open("Organism_group/node_list_temp/"+ organism_name[i] +"_nodelist", "r")
        f_grp_open = open("Organism_group/group_based_on_identity_score/" + organism_name[i] +"_Group", "w")

        G = nx.Graph()

        ######################### Reading nodelist input file ###
        while(True):
            line = f_open_nodelist.readline()
            if(line[0:18] == "No of input files:"):
                break

            G.add_node(line.strip())
        
        line = line.split(" ")   
        no_fo_input_files = int(line[4].strip())


        for file_no in range(1, no_fo_input_files+1):
            f_open = open("Organism_identity_score/identity_score_list_temp/"+ organism_name[i] + "_RNA_chain_list_" + str(file_no) + "_Identity_Score_list", "r")
            f_identity_grp_open = open("Organism_group/identity_100_list/" + organism_name[i] +"_identity_100", "w")

            
            line = f_open.readline()
            while(True):
                line = f_open.readline()
                if(line == ''):
                    break
                line = line.strip().split("\t")
                pdb1 = line[0]
                chain1 = line[1]
                pdb2 = line[3]
                chain2 = line[4]
                identity = line[6]

                if(float(identity)) >= identity_threshold:
                    node1 = pdb1 + "_" + chain1
                    node2 = pdb2 + "_" + chain2
                    G.add_edge(node1, node2)
                    if(float(identity)) == 100:
                        f_identity_grp_open.write(node1 + "\t" + node2 + "\n")
                    
                    
            f_open.close()

        connected_graph_list = list(nx.connected_components(G))

        for component in range(0,len(connected_graph_list)):
            f_grp_open.write("Group no: "+str(component) + "\n")
            f_grp_open.write(str(connected_graph_list[component])+ "\n")
            
            

        f_grp_open.close()
        f_open_nodelist.close()

