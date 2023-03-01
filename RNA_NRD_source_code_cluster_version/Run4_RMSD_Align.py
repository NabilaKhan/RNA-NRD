import math
import os
import shutil
import sys
import argparse
import glob
import logging
import time




def generate_RMSD_align(max_chain):

    print("\n ------- Generate RNA 3D strcture alignment ------- \n")

    os.chdir('src/')
    
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

        # Count no of clusters for each organism to understand number of batch files for each organism
        f_open = open("Organism_group/group_based_on_identity_score/" + organism_name[i] + "_Group", "r")
        no_of_input_files = 0
        cluster_length_list = []

        while(True):
            
            line = f_open.readline()
            if(line == ''):
                break
            line = f_open.readline()
            line = line.strip("\n").lstrip("{").rstrip("}")
            line_list = line.split(',')
            cluster_length_list.append(len(line_list))
            no_of_input_files += 1


        
        batch_info_dir = "Organism_group/batch_file_list_temp/"
        os.chdir(batch_info_dir)
        
       
        for input_file_ind in range(0, no_of_input_files):

            no_of_chains = cluster_length_list[input_file_ind]
            list_no = 1

            while(no_of_chains > 0):
              
                batch_file_name =  "Align_" + organism_name[i] + "_batch_file_cluster_" + str(input_file_ind+1) + "_list_" + str(list_no)
                os.system("chmod +x " + batch_file_name)
                os.system("sbatch " + batch_file_name)
                list_no += 1
                no_of_chains = no_of_chains - max_chain
                
            if(no_of_chains > 0):
        
                batch_file_name =  "Align_" + organism_name[i] + "_batch_file_cluster_" + str(input_file_ind+1) + "_list_" + str(list_no)
                os.system("chmod +x " + batch_file_name)
                os.system("sbatch " + batch_file_name)

        os.chdir("../..")
            
        f_open.close()
    

def main():

    process_start_time = time.time()
    parser = argparse.ArgumentParser(description='Generate alignment for RNA 3D structures')
    parser.add_argument('-fr', nargs='?', default=50, const=50, help="Number of RNA chains per file for calculating RMSD. It's used to divide total number of RNA chains in an organism into smaller sections to make parallel processing faster while calculating RMSD. Default: 50.")

    try:
        args = parser.parse_args()
    except Exception as e:
        parser.print_help()
        sys.exit()


    rmsd_rna_chains_per_file = int(args.fr)
    

    generate_RMSD_align(rmsd_rna_chains_per_file)


if __name__ == '__main__':
    main()
