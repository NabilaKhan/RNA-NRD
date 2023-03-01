import math
import os
import shutil
import sys
sys.path.append('src/')

  
def generate_file_for_RMSD_calculation(RNA_chains_per_file):

    print("\n ------- Generate necessary files and sbatch files for aligning RNA 3D strctures ------- \n")
    
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
     

    ######################### Creating Temporary Folders for Calculating RMSD from Complete Graph Parallely #########################
    path = "Organism_group/" + "input_file_list_temp"

    if(os.path.isdir(path)):
        shutil.rmtree(path)    
    os.mkdir(path)


    path_batch_file = "Organism_group/" + "batch_file_list_temp"

    if(os.path.isdir(path_batch_file)):
        shutil.rmtree(path_batch_file)    
    os.mkdir(path_batch_file)


    path_rmsd = "Organism_RMSD/" + "rmsd_cluster_tmp"

    if(os.path.isdir(path_rmsd)):
        shutil.rmtree(path_rmsd)    
    os.mkdir(path_rmsd)


    path_rmsd_cluster = "Organism_group/" + "cluster_based_on_RMSD_temp"

    if(os.path.isdir(path_rmsd_cluster)):
        shutil.rmtree(path_rmsd_cluster)    
    os.mkdir(path_rmsd_cluster)


    path_log = "Organism_group/" + "Log_file"

    if(os.path.isdir(path_log)):
        shutil.rmtree(path_log)    
    os.mkdir(path_log)
    ######################### Done Creating Temporary Folders #########################



    for i in range(0, len(organism_name)):

        f_open = open("Organism_group/group_based_on_identity_score/" + organism_name[i] + "_Group", "r")
        
        chain_list = {}
        key_val = 0
        list_count = 0

        ######################### For each organism, write the nodes in each cluster to different files
        while(True):
            
            line = f_open.readline()
            if(line == ''):
                break
            line = f_open.readline()
            line = line.strip("\n").lstrip("{").rstrip("}")
            line_list = line.split(',')
            list_count += 1

            f_open_input_file = open("Organism_group/input_file_list_temp/" + organism_name[i] + "_cluster_" + str(list_count), "w")
            

            for ch in range(0, len(line_list)):
                line_list[ch] = line_list[ch].strip().lstrip("'").rstrip("'")
                f_open_input_file.write(line_list[ch] + "\n")   

            f_open_input_file.close()


            len_cluster_list = len(line_list)
            if(len_cluster_list % RNA_chains_per_file) == 0:
                no_of_cluster_inputs = math.floor(len_cluster_list / RNA_chains_per_file)
            else:
                no_of_cluster_inputs = math.floor(len_cluster_list / RNA_chains_per_file) + 1

          
            
            cluster_list_count = 0
            for fl in range(0, no_of_cluster_inputs):

                f_open_input_file_list = open("Organism_group/input_file_list_temp/" + organism_name[i] + "_cluster_" + str(list_count) + "_list_" + str(fl+1), "w")
                for list_ch in range(fl * RNA_chains_per_file, (fl+1) * RNA_chains_per_file):
                    if(list_ch == len_cluster_list):
                        break
                    f_open_input_file_list.write(line_list[list_ch] + "\n") 
                f_open_input_file_list.close()
                

                ######################### Write Batch File for preprocessing in generation of RMSD from clusters/gropus of chains ###
                f_batch_file_pre = open("Organism_group/batch_file_list_temp/Preprocess_"+ organism_name[i] +"_batch_file_cluster_" + str(list_count) + "_list_" + str(fl+1) , "w")
                f_batch_file_pre.write("#!/bin/bash\n")
                f_batch_file_pre.write("#SBATCH --output=Preprocess_" + organism_name[i] + "_batch_file_output_" + str(list_count) + "_list_" + str(fl+1) + "\n")
                f_batch_file_pre.write("#SBATCH --cpus-per-task=1\n")
                f_batch_file_pre.write("#SBATCH --mem=16G\n\n\n")

                new_cluster_file_name = organism_name[i] +"_cluster_" + str(list_count)
                new_input_file_name = organism_name[i] + "_cluster_" + str(list_count) + "_list_" + str(fl+1)
                f_batch_file_pre.write("python3 ../../Generate_RMSD_Preprocess_from_Connected_Graph_Parallel_Batch_file.py " + organism_name[i] + " " + new_cluster_file_name + " " + new_input_file_name + " " + str(RNA_chains_per_file) + " " + str(fl+1) + " \n")
                f_batch_file_pre.close()

                
                ######################### Write Batch File for alignment in generation of RMSD from clusters/gropus of chains ###
                f_batch_file = open("Organism_group/batch_file_list_temp/Align_"+ organism_name[i] +"_batch_file_cluster_" + str(list_count) + "_list_" + str(fl+1) , "w")
                f_batch_file.write("#!/bin/bash\n")
                f_batch_file.write("#SBATCH --output=Align_" + organism_name[i] + "_batch_file_output_" + str(list_count) + "_list_" + str(fl+1) + "\n")
                f_batch_file.write("#SBATCH --cpus-per-task=1\n")
                f_batch_file.write("#SBATCH --mem=16G\n\n\n")

                new_cluster_file_name = organism_name[i] +"_cluster_" + str(list_count)
                new_input_file_name = organism_name[i] + "_cluster_" + str(list_count) + "_list_" + str(fl+1)
                f_batch_file.write("python3 ../../Generate_RMSD_Align_from_Connected_Graph_Parallel_Batch_file.py " + organism_name[i] + " " + new_cluster_file_name + " " + new_input_file_name + " " + str(RNA_chains_per_file) + " " + str(fl+1) + " \n")
                f_batch_file.close()


        f_open.close()
