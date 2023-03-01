import argparse
import sys
sys.path.append('src/')
import os
import glob
import logging
import time


from Cluster_Based_on_RMSD_Value import *
from Selecting_Representative import *



def main():

    print(" ------- Generate Non-redudant dataset along with representative for each cluster (with organism division) ------- \n")
    
    os.chdir('src/')

    process_start_time = time.time()
    parser = argparse.ArgumentParser(description='Generate Non-redundant Database')
    parser.add_argument('-o', nargs='?', default='Output/RNA-NRD_dataset', const='Output/RNA-NRD_dataset', help="Output file name and location. Default: 'Output/RNA-NRD_dataset'.")
    parser.add_argument('-r', nargs='?', default=4, const=4, help='Provide RMSD threshold. Default: 4.')
    parser.add_argument('-a', nargs='?', default=80, const=80, help='Provide alignment ratio percent threshold. Default: 80.')
    parser.add_argument('-fr', nargs='?', default=50, const=50, help="Number of RNA chains per file for calculating RMSD. It's used to divide total number of RNA chains in an organism into smaller sections to make parallel processing faster while calculating RMSD. Default: 50.")

    try:
        args = parser.parse_args()
    except Exception as e:
        parser.print_help()
        sys.exit()


    output_dir = args.o
    rmsd_threshold = float(args.r)
    align_ratio_threshold = int(args.a)
    rmsd_rna_chains_per_file = int(args.fr)
    

    generate_cluster_based_on_RMSD(rmsd_threshold, align_ratio_threshold, rmsd_rna_chains_per_file)
    generate_output(output_dir, rmsd_rna_chains_per_file)


if __name__ == '__main__':
    main()
