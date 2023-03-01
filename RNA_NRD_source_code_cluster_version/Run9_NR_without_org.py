import argparse
import sys
sys.path.append('src/')
import os
import glob
import logging
import time


from RNA_NRD_without_org_div import *




def main():

    print("\n ------- Generate Non-redudant dataset along with representative for each cluster (without organism division) ------- \n")
    
    os.chdir('src/')

    process_start_time = time.time()
    parser = argparse.ArgumentParser(description='Generate RNA-NRD-without-Organism-Division Dataset')
    parser.add_argument('-o', nargs='?', default='RNA-NRD_without_Organism_Division_dataset', const='RNA-NRD_without_Organism_Division_dataset', help="Output file name and location. Default: 'RNA-NRD_without_Organism_Division_dataset'.")
    
    try:
        args = parser.parse_args()
    except Exception as e:
        parser.print_help()
        sys.exit()


    output_dir = args.o

    generate_output_without_org(output_dir)


if __name__ == '__main__':
    main()
