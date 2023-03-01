import argparse
import sys
sys.path.append('src/')
import os
import glob
import logging
import time


from Merge_Organism_SI import *
from Merge_Organism_RMSD import *
from Merge_Organism_Members import *
from RNA_NRD_without_org_div import *




def main():

    print("\n ------- Removing organism based division ------- \n")

    print("\n ------- Generating sequence alignement among cluster representatives ------- \n")
    
    os.chdir('src/')

    process_start_time = time.time()
    parser = argparse.ArgumentParser(description='Generate RNA-NRD-without-Organism-Division Dataset')
    parser.add_argument('-s', nargs='?', default=80, const=80, help='Provide sequence identity threshold. Default: 80.')
    parser.add_argument('-r', nargs='?', default=4, const=4, help='Provide RMSD threshold. Default: 4.')
    parser.add_argument('-a', nargs='?', default=80, const=80, help='Provide alignment ratio percent threshold. Default: 80.')
    parser.add_argument('-i', nargs='?', default='RNA-NRD_dataset', const='RNA-NRD_dataset', help="RNA-NRD dataset name and location. Default: 'RNA-NRD_dataset'.")


    
    try:
        args = parser.parse_args()
    except Exception as e:
        parser.print_help()
        sys.exit()


    rmsd_threshold = float(args.r)
    align_ratio_threshold = int(args.a)
    seq_identity_threshold = int(args.s)
    nr_dataset = args.i
    

    align_RP_seq(seq_identity_threshold, nr_dataset)



if __name__ == '__main__':
    main()
