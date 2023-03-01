import argparse
import sys
sys.path.append('src/')
import os
import glob
import logging
import time

from Select_RNA_chains import *
from Organism_specific_RNA_chains import *
from Generate_Parallel_Batch_File_for_Identity_Score import *



def main():

    process_start_time = time.time()
    parser = argparse.ArgumentParser(description='Preprocess RNA chains')
    parser.add_argument('-i', nargs='?', default='Data/Data.txt', const='Data/Data.txt', help="Input file location containing RNA chains. Default: 'Data/Data.txt'.")
    parser.add_argument('-fi', nargs='?', default=100, const=100, help="Number of RNA chains per file for calculating sequence identity. It's used to divide total number of RNA chains in an organism into smaller sections to make parallel processing faster while calculating sequence identity. Default: 100")
    
    try:
        args = parser.parse_args()
    except Exception as e:
        parser.print_help()
        sys.exit()

    user_input_fname = args.i
    sequence_identity_rna_chains_per_file = int(args.fi)  ### Number of RNA chains in each batch file

    select_RNA_chains(user_input_fname)
    divide_chains_based_on_organism()
    generate_file_for_identity_score_calculation(sequence_identity_rna_chains_per_file)



if __name__ == '__main__':
    main()
