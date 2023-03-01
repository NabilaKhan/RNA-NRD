import math
import shutil
import os
import argparse
import glob
import logging
import time
import sys
##sys.path.append('src/')


def generate_identity_score(RNA_chains_per_file):

    print("\n ------- Generating sequence identity score ------- \n")

    os.chdir('src/')
    #print(os.getcwd())

    f_organism_name = open("Organism_list/Current_Organism_list","r")
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

        
        len_chain_list = len(chain_list)
        no_of_input_files = math.ceil(len_chain_list/RNA_chains_per_file)

        batch_info_dir = "Organism_identity_score/batch_file_list_temp/"
        os.chdir(batch_info_dir)
        
       
        for input_file_ind in range(0, no_of_input_files):
            
            batch_file_name =  organism_name[i] + "_batch_file_" + str(input_file_ind+1)
            os.system("chmod +x " + batch_file_name)
            os.system("sbatch " + batch_file_name)

        os.chdir("../..")
            
        f_open.close()



def main():

    process_start_time = time.time()
    parser = argparse.ArgumentParser(description='Generate Sequence Identity Score')
    parser.add_argument('-fi', nargs='?', default=100, const=100, help="Number of RNA chains per file for calculating sequence identity. It's used to divide total number of RNA chains in an organism into smaller sections to make parallel processing faster while calculating sequence identity.")

    try:
        args = parser.parse_args()
    except Exception as e:
        parser.print_help()
        sys.exit()

    sequence_identity_rna_chains_per_file = int(args.fi)  ### Number of RNA chains in each batch file


    generate_identity_score(sequence_identity_rna_chains_per_file)



if __name__ == '__main__':
    main()

    


