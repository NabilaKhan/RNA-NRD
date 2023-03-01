import subprocess
import os.path
import sys
sys.path.append('src/')
import glob
import time


organism_name = sys.argv[1]
organism_cluster_file_name = sys.argv[2]
organism_input_file_name = sys.argv[3]
RNA_chains_per_file = int(sys.argv[4])
organism_file_no = int(sys.argv[5])


organism_cluster_file_path = 'Organism_group/input_file_list_temp/' + organism_cluster_file_name
organism_input_file_path = 'Organism_group/input_file_list_temp/' + organism_input_file_name
print(organism_name)
print(organism_input_file_name)


f_all_chains_in_cluster = open('../../' + organism_cluster_file_path, "r")
f_current_cluster = open('../../' + organism_input_file_path, "r")
f_out = open('../../' + 'Organism_RMSD/rmsd_cluster_tmp/' + organism_input_file_name, "w")
log_file = open('../../' + 'Organism_group/Log_file/Preprocess_' + organism_input_file_name, "w")
start_time = time.time()


######################### List of all chains in the cluster of an organism
Total_chain_list = {}
tot_key_val = 0
while(True):
    
    line = f_all_chains_in_cluster.readline()
    if(line == ''):
        break
    line = line.split('_')
    pdb_id, chain,  = line[0], line[1].strip()
    Total_chain_list[tot_key_val] = [pdb_id, chain]
    
    tot_key_val += 1

len_tot_chain_list = len(Total_chain_list)


######################### List of chains currently considered of an organism in the cluster
comparing_chain_list = {}
key_val = 0

while(True):
    
    line = f_current_cluster.readline()
    if(line == ''):
        break
    line = line.split('_')
    pdb_id, chain,  = line[0], line[1].strip()
    comparing_chain_list[key_val] = [pdb_id, chain]
    key_val += 1

len_chain_list = len(comparing_chain_list)
comparing_key_start = (organism_file_no - 1) * RNA_chains_per_file


######################### Preprocessing pdb/cif files using STAR3D_dssr
for main_key in comparing_chain_list:
    
    comparing_chain = comparing_chain_list[main_key]
    pdb_id = comparing_chain[0]         
    chain = comparing_chain[1]
    print("Preprocessing " + pdb_id + ' ' + chain)

    STAR3D_path = '../' + '../STAR3D_source_dssr/'
    p1 = subprocess.call('java -cp STAR3D.jar Preprocess ' + pdb_id + ' ' + chain, shell = True, cwd = STAR3D_path)

    

f_current_cluster.close()
f_out.close()
run_time = round((time.time() - start_time)/3600,2)
log_file.write("--- %s hours ---" % (str(run_time)))
