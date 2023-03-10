from Bio.Seq import Seq 
##from Bio.pairwise2 import format_alignment
##from Bio import pairwise2
from Bio.Align import PairwiseAligner
import sys
sys.path.append('src/')


organism_name = sys.argv[1]
organism_input_file_name = sys.argv[2]
organism_file_no = int(sys.argv[3])
RNA_chains_per_file = int(sys.argv[4])


organism_input_file_path = 'Organism_identity_score/input_file_list_temp/' + organism_input_file_name


f_open = open('../../' + "Organism_chains/" + organism_name + "_RNA_chains", "r")
f_current_chain_open = open('../../' + organism_input_file_path, "r")
f_out = open('../../' + "Organism_identity_score/identity_score_list_temp/" + organism_input_file_name + "_Identity_Score_list", "w")
f_identity = open('../../' + "Organism_group/identity_100_list/"+ organism_name + "_identity_100", "w")


######################### List of all chains of an organism
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

    

######################### List of chains currently considered of an organism
comparing_chain_list = {}
line = f_current_chain_open.readline()
key_val = 0

while(True):
    
    line = f_current_chain_open.readline()
    if(line == ''):
        break
    line = line.split('\t')
    pdb_id, chain, RNA_name, sequence = line[0], line[1].strip(), line[2], line[3]
    chain = chain.replace("\"", "")
    comparing_chain_list[key_val] = [pdb_id, chain, RNA_name, sequence]
    key_val += 1


len_chain_list = len(chain_list)
comparing_key_start = (organism_file_no - 1) * RNA_chains_per_file

f_out.write("PDB_id1\tchain1\tRNA_name1\tPDB_id2\tchain2\tRNA_name2\tIdentity Score\n")


  
for key in comparing_chain_list:

    current_chain = comparing_chain_list[key]
    seq1 = current_chain[3].strip()
    seq1_len = len(seq1)
    pdb_id1 = current_chain[0]
    chain1 = current_chain[1]
    RNA_name1 = current_chain[2]
    comparing_key = comparing_key_start + key

    chain1 = chain1.replace(" ", "")
    chain1_list = chain1.split(",")

    fix_identity = 100

    for ch1 in range(0, len(chain1_list)):
        for c in range(ch1+1, len(chain1_list)):
            f_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%0.3f\n" % (pdb_id1, chain1_list[ch1], RNA_name1, pdb_id1, chain1_list[c], RNA_name1, fix_identity))
            f_identity.write(pdb_id1 + '_' + chain1_list[ch1] + '\t' + pdb_id1 + '_' + chain1_list[c] + '\n')
    
    for key2 in range(comparing_key + 1, len_chain_list):

        comparing_chain = chain_list[key2]
        pdb_id2 = comparing_chain[0]
        seq2 = comparing_chain[3].strip()
        seq2_len = len(seq2)
        chain2 = comparing_chain[1]
        RNA_name2 = comparing_chain[2]

        chain2 = chain2.replace(" ", "")
        chain2_list = chain2.split(",")

          
        if((seq1_len < 2*seq2_len) and (seq2_len < 2*seq1_len)):
         
            identity = 0.0
##            alignments = pairwise2.align.globalxx(seq1, seq2)
            aligner = PairwiseAligner(mode='global', match_score=1, mismatch_score=0)
            alignments = aligner.align(seq1, seq2)
            
            for alignment in alignments:
                
##                score = str(alignment).split()
##                score = float(score[2].strip(','))
                score = float(alignment.score)
                                
                if(seq1_len >= seq2_len):
                    identity = (score/seq1_len) * 100
                else:
                    identity = (score/seq2_len) * 100

                for ch1 in range(0, len(chain1_list)):
                    for ch2 in range(0, len(chain2_list)):
                        f_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%0.3f\n" % (pdb_id1, chain1_list[ch1], RNA_name1, pdb_id2, chain2_list[ch2], RNA_name2, identity))
                        if(identity == 100):
                            f_identity.write(pdb_id1 + '_' + chain1_list[ch1] + '\t' + pdb_id2 + '_' + chain2_list[ch2] + '\n')
                break     

f_open.close()
f_current_chain_open.close()
f_out.close()
