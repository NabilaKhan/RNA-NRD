import math
import os
import subprocess

pdbx_dir = 'input_files/'
pdb_fasta_mapping_dir = 'pdb_fasta/'
#pdb_chains = {'6ypu': ['2','4'], '6ytf': ['3']}
#pdb_chains = {'6ypu': ['2','4']}
#pdb_chains = {'6ypu': ['2'], '6ytf': ['3']}


#organism_division = ['ABC', 'DEFGH', 'IJKLM', 'NOPQR', 'STUV', 'WXYZ']
#organism_division = ['ABC']
#organism_division = ['DEFGH']
#organism_division = ['IJKLM']
#organism_division = ['NOPQR']
#organism_division = ['STUV']
organism_division = ['WXYZ']



f_pairwise_rmsd = open(organism_division[0] + "-rmsd","r")
f_20_rmsd = open(organism_division[0] + "_20_rmsd","w")
count_20  = 0

while(True):
    line = f_pairwise_rmsd.readline()
    if(line == ''):
        break

    line = line.strip().split("\t")
    pdb1, chain1, pdb2, chain2, rmsd = line[0], line[1], line[2], line[3], line[4]
    if(float(rmsd) == 20):
        print(pdb1, chain1, pdb2, chain2, rmsd)
        f_20_rmsd.write("%s\t%s\t%s\t%s\t%s\n" % (pdb1, chain1, pdb2, chain2, rmsd))
        count_20+= 1

print(count_20)


f_pairwise_rmsd.close()
f_20_rmsd.close()






   


