import argparse
import sys
sys.path.append('src/')
import os
import glob
import logging
import time

##
##from Merge_Organism_SI import *
##from Merge_Organism_RMSD import *
##from Merge_Organism_Members import *
##from RNA_NRD_without_org_div import *




def main():

    print("\n ------- Generating structure alignement among cluster representatives ------- \n")
    
    #os.chdir('src/')

    process_start_time = time.time()
    parser = argparse.ArgumentParser(description='Generate RNA-NRD-without-Organism-Division Dataset')
    parser.add_argument('-r', nargs='?', default=4, const=4, help='Provide RMSD threshold. Default: 4.')
    parser.add_argument('-a', nargs='?', default=80, const=80, help='Provide alignment ratio percent threshold. Default: 80.')
    parser.add_argument('-i', nargs='?', default='RNA-NRD_dataset', const='RNA-NRD_dataset', help="RNA-NRD dataset name and location. Default: 'RNA-NRD_dataset'.")


    
    try:
        args = parser.parse_args()
    except Exception as e:
        parser.print_help()
        sys.exit()


    rmsd_threshold = float(args.r)
    align_ratio_threshold = float(args.a)
    nr_dataset = args.i


    ######################### Write Batch File for RP RMSD #########################
    batch_file_name = "src/Organism_merge/Batch_file/RP_RMSD_batch"
    f_batch_file = open(batch_file_name , "w")
    f_batch_file.write("#!/bin/bash\n")
    f_batch_file.write("#SBATCH --output=RP_RMSD_batch_file_output\n")
    f_batch_file.write("#SBATCH --cpus-per-task=1\n")
    f_batch_file.write("#SBATCH --mem=16G\n\n\n")

    f_batch_file.write("python3 ../../Merge_Organism_RMSD.py " + str(rmsd_threshold) + " " + str(align_ratio_threshold) + " " + nr_dataset + " \n")
    f_batch_file.close()


    ######################### Run Batch File for RP RMSD #########################
    batch_info_dir = "src/Organism_merge/Batch_file/"
    os.chdir(batch_info_dir)
    #os.system("chmod +x " + batch_file_name)
    os.system("sbatch RP_RMSD_batch")
        
    



if __name__ == '__main__':
    main()
