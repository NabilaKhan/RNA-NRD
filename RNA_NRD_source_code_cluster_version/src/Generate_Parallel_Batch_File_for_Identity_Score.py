import math
import os
import shutil
import sys
sys.path.append('src/')


def generate_file_for_identity_score_calculation(RNA_chains_per_file):

    print("\n ------- Generating necessary files and sbatch files for calculating RNA sequence alignment ------- \n")

    
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

    

    ######################### Creating Temporary Folders for Calculating Identity Score Paralelly #########################
    path = "Organism_identity_score/" + "input_file_list_temp"

    if(os.path.isdir(path)):
        shutil.rmtree(path)    
    os.mkdir(path)


    path_batch_file = "Organism_identity_score/batch_file_list_temp"

    if(os.path.isdir(path_batch_file)):
        shutil.rmtree(path_batch_file)    
    os.mkdir(path_batch_file)


    path_identity = "Organism_identity_score/identity_score_list_temp"

    if(os.path.isdir(path_identity)):
        shutil.rmtree(path_identity)    
    os.mkdir(path_identity)


    path_nodelist = "Organism_group/node_list_temp"

    if(os.path.isdir(path_nodelist)):
        shutil.rmtree(path_nodelist)    
    os.mkdir(path_nodelist)


    path_nodelist = "Organism_group/group_based_on_identity_score"

    if(os.path.isdir(path_nodelist)):
        shutil.rmtree(path_nodelist)    
    os.mkdir(path_nodelist)
    ######################### Done Creating Temporary Folders #########################



    ######################### Creating Temporary Folders for RNA chains with 100 identiy score #########################
    path_100_identity = "Organism_group/" + "identity_100_list"

    if(os.path.isdir(path_100_identity)):
        shutil.rmtree(path_100_identity)    
    os.mkdir(path_100_identity)
    ######################### Done Creating Temporary Folders #########################



    for i in range(0, len(organism_name)):

        f_open = open("Organism_chains/" + organism_name[i] + "_RNA_chains", "r")
        
        chain_list = {}
        line = f_open.readline()
        key_val = 0

        while(True):
            
            line = f_open.readline()
            if(line == ''):
                break
            line = line.split('\t')
            
            pdb_id, chain, RNA_name, sequence = line[0], line[1].strip(), line[9], line[13]
            chain = chain.replace("\"", "")
            chain_list[key_val] = [pdb_id, chain, RNA_name, sequence]
            key_val += 1

        f_open.close()
        

        len_chain_list = len(chain_list)
        no_of_input_files = math.ceil(len_chain_list/RNA_chains_per_file)
        

        ######################### Writing List of Nodes to temp folder 'node_list_temp' #########################
        f_nodelist_out = open("Organism_group/node_list_temp/" + organism_name[i] + "_nodelist", "w")
        
        for key in chain_list:
            chain_id = chain_list[key][1]
            chain_id = chain_id.replace(" ", "")
            chain_id_list = chain_id.split(",")

            for ch in range(0, len(chain_id_list)):
                node = chain_list[key][0] + "_" + chain_id_list[ch]
                f_nodelist_out.write(node+"\n")
                
        f_nodelist_out.write("No of input files: "+ str(no_of_input_files) +"\n")

        f_nodelist_out.close()
        ######################### Done Writing #########################
        
            
        for input_file_ind in range(0, no_of_input_files):
            
            ######################### Write parallel input file list in the identity_score folder #########################
            f_out = open("Organism_identity_score/input_file_list_temp/"+ organism_name[i] +"_RNA_chain_list_" + str(input_file_ind+1) , "w")
            f_out.write("PDB_ID\tSequence\tMacromolecule_Name\tChain_ID\n")
            start_key_index = input_file_ind * RNA_chains_per_file + 0
            if(((input_file_ind + 1) * RNA_chains_per_file) <= len_chain_list):
                end_key_index = start_key_index + RNA_chains_per_file
            else:
                end_key_index = start_key_index + len_chain_list % RNA_chains_per_file


            for key_index in range(start_key_index, end_key_index):
                
                current_chain = chain_list[key_index]
                PDB_id = current_chain[0]
                Chain_id = current_chain[1]
                Macromolecule_name = current_chain[2]
                Sequence = current_chain[3]
                f_out.write("%s\t%s\t%s\t%s" % (PDB_id, Chain_id, Macromolecule_name, Sequence))
            
            f_out.close()


            ######################### Write Batch File for identity score #########################
            f_batch_file = open("Organism_identity_score/batch_file_list_temp/"+ organism_name[i] +"_batch_file_" + str(input_file_ind+1) , "w")
            f_batch_file.write("#!/bin/bash\n")
            f_batch_file.write("#SBATCH --output=" + organism_name[i] + "_batch_file_output_" + str(input_file_ind+1) + "\n")
            f_batch_file.write("#SBATCH --cpus-per-task=1\n")
            f_batch_file.write("#SBATCH --mem=16G\n\n\n")

            new_RNA_list_file_name = organism_name[i] +"_RNA_chain_list_" + str(input_file_ind+1)
            f_batch_file.write("python3 ../../Generate_Identity_Score.py " + organism_name[i] + " " + new_RNA_list_file_name + " " + str(input_file_ind+1) + " " + str(RNA_chains_per_file) + " \n")
            f_batch_file.close()

        f_open.close()
            


