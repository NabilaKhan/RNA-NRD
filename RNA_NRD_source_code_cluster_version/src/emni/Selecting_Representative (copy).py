import math
import os
import shutil
import networkx as nx
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import sys
import time



def generate_output(output_dir, max_chain):

    print(" ------- Generate Non-redudant dataset along with representative for each cluster ------- \n")


    group_id = 1
    now = time.gmtime()
    current_time = time.mktime(now)

    final_table = open("../" + output_dir,"w")
    final_table.write("Cluster_ID\tRepresentative\tRedundant_cluster\tOrganism\tMacromolecule_name\tFamily_name\n")
    
    f_organism_name = open("Organism_list/Current_Organism_list", "r")
    organism_name = []
    organism_count = 0

    while(True):
        line = f_organism_name.readline()
        
        if line == "":
            break
        
        line = line.strip("\n")
        organism = line.replace(" ", "_")
        organism_name.append(organism)
        organism_count += 1

    f_organism_name.close()  



    for i in range(0, len(organism_name)):

        f_chain = open("Organism_chains/" + organism_name[i] + "_RNA_chains", "r")
        chain_info = {}

        while(True):
            chain_line = f_chain.readline()
            if(chain_line == ""):
                break

            chain_line = chain_line.split("\t")
            pdb, chain, chain_len, resolution, release_date, exp_method, macromolecule_name, fam_name, num_coord, num_bp = chain_line[0], chain_line[1], chain_line[4], chain_line[5], chain_line[7], chain_line[8], chain_line[9], chain_line[10], chain_line[11], chain_line[12]
            info_list = [chain_len, resolution, release_date, exp_method, macromolecule_name, fam_name, num_coord, num_bp]
            pdb_chain = pdb + '_' + chain
            chain_info[pdb_chain] = info_list

      
        ######################### Count no of clusters for each organism to understand number of batch files for each organism
        f_open = open("Organism_group/group_based_on_identity_score/" + organism_name[i] + "_Group", "r")
        no_of_input_files = 0
        cluster_length_list = []

        while(True):
            
            line = f_open.readline()
            if(line == ''):
                break
            line = f_open.readline()
            line = line.strip("\n").lstrip("{").rstrip("}")
            line = line.replace(" ", "").replace("'","")
            line_list = line.split(',')
            cluster_length_list.append(len(line_list))
            no_of_input_files += 1
        
       
        for input_file_ind in range(0, no_of_input_files):
                
            f_open_group = open("Organism_group/cluster_based_on_RMSD_temp/" + organism_name[i] + "_cluster_" + str(input_file_ind+1), "r")

           
            while(True):
                line_grp = f_open_group.readline()
                if(line_grp == ''):
                    break

                line_grp = f_open_group.readline()
                line_grp = line_grp.strip("\n").lstrip("{").rstrip("}").replace(" ", "").replace("'","")
                line_list = line_grp.split(',')
              
                line_grp = f_open_group.readline()
                line_grp = line_grp.strip("\n").lstrip("[").rstrip("]").replace(" ", "")
                degree_list = line_grp.split(',')
                
                macromolecule_list = []
                family_list = []
                release_list = []
                coord_list = []
                bp_list = []
                resolution_list = []
                chain_len_list = []
                exp_list = []

                quality_score_list = []
                
                for li in range(0,len(line_list)):
                    
                    macromolecule_list.append(chain_info[line_list[li]][4])
                    fam_temp = chain_info[line_list[li]][5].lstrip("[").rstrip("]")
                    fam_temp = fam_temp.strip("'")
                    family_list.append(fam_temp)
                    

                    if chain_info[line_list[li]][1] == '-':
                        resolution = 0
                    elif '_' in chain_info[line_list[li]][1]:
                        resolution = chain_info[line_list[li]][1]
                        resolution = resolution.split('_')[0].strip(',')
                    else:
                        resolution = chain_info[line_list[li]][1]

                    temp = (int(chain_info[line_list[li]][6]) + int(chain_info[line_list[li]][7]))/100 + (float(resolution)/10) + (int(degree_list[li])/20)
                    quality_score_list.append(temp)

                ######################### Generate Representative with maximum quality score
                max_value = max(quality_score_list)
                max_index_list = []
                for k in range(0, len(quality_score_list)):
                    if quality_score_list[k] == max_value:
                        max_index_list.append(k)
               
                ######################### Tie breaker
                ######################### Tie based on quality Score
                if len(max_index_list) > 1:

                    ### Tie breaker based on release date
                    release_list = []
                    for ind in range(0, len(max_index_list)):
                        chain_date = chain_info[line_list[max_index_list[ind]]][2]
                        date_form = time.strptime(chain_date, "%Y-%m-%d")
                        chain_time = time.mktime(date_form)
                        difference = (current_time - chain_time) / 3600
                        release_list.append(difference)


                    ### Generate Representative with latest date
                    min_date = min(release_list)
                    min_date_index_list = []
                    for k in range(0, len(release_list)):
                        if release_list[k] == min_date:
                            min_date_index_list.append(max_index_list[k])

                  
                    ### Still tie on release date 
                    if len(min_date_index_list) > 1:

                        ### Tie breaker based on exp method
                        exp_list = []
                        for ex in range(0, len(min_date_index_list)):
                            exp_list.append(chain_info[line_list[min_date_index_list[ex]]][3])
                            
                        ### Select Representative with X-RAY_DIFFRACTION method if there is any
                        exp_index_list = []
                        if 'X-RAY_DIFFRACTION' in exp_list:
                            for exi in range(0, len(exp_list)):
                                if exp_list[exi] == 'X-RAY_DIFFRACTION':
                                    exp_index_list.append(min_date_index_list[exi])
                        else:
                            exp_index_list = min_date_index_list

                       
                        ### Still tie on experiment method 
                        if len(exp_index_list) > 1:

                            ### Tie breaker based on chain length
                            len_list = []
                            for l in range(0, len(exp_index_list)):
                                len_list.append(chain_info[line_list[exp_index_list[l]]][0])
                               
                            ### Select Representative based on max chain length
                            max_len = max(len_list)
                            max_chain_index_list = []
                            for ch in range(0, len(len_list)):
                                if len_list[ch] == max_len:
                                    max_chain_index_list.append(exp_index_list[ch])
                                    
                            ### If Still tie on chain length, select the first one from the multiple ones with the longest length; otherwise
                            ### selecting representative based on minimum quality score, latest release date, experiment method and longer chain length    
                            representative_ind = max_chain_index_list[0]    

                        ### Selected representative based on minimum quality score, latest release date and experiment method        
                        else:
                            representative_ind = exp_index_list[0]                                 
                    ### Selected representative based on minimum quality score and latest release date
                    else:
                        representative_ind = min_date_index_list[0]
                        
                ### Selected representative based on minimum quality score
                else:
                    representative_ind = max_index_list[0]

            
                final_table.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (group_id, line_list[representative_ind], ','.join(line_list), organism_name[i], ','.join(list(set(macromolecule_list))), ','.join(list(set(family_list)))))
                group_id += 1

            f_open_group.close() 
        f_open.close()
    final_table.close()



