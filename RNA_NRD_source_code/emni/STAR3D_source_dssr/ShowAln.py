'''
This script shows the 3D align result in PyMol.
Make sure PyMol is installed in your system.
usage: 
1. open PyMol
2. import ShowAln into PyMol: run /path/to/script/ShowAln.py
3. show the alignment: ShowAln aln_pdb, aln_file
	aln_pdb: the output aligned pdb file
	aln_file: the output aligned file 
'''

from pymol import cmd, stored

def ShowAln(aln_pdb, aln_file):
	aln_index1=[]
	aln_index2=[]
	fp=open(aln_file)
	for line in fp.readlines():
		if line.startswith('#'):
			continue
		line=line.strip().split('<->')
		aln_index1.append(line[0][2:])
		aln_index2.append(line[1][2:])
	fp.close()

	cmd.load(aln_pdb, "star3d_aln")
	cmd.split_states("star3d_aln")
	cmd.set_name("star3d_aln_0001", "aln_seq1")
	cmd.set_name("star3d_aln_0002", "aln_seq2")
	
	cmd.show("ribbon", 'aln_seq1')
	cmd.hide('lines', 'aln_seq1')
	cmd.hide('spheres', 'aln_seq1')
	cmd.hide('sticks', 'aln_seq1')
	
	cmd.show("ribbon", "aln_seq2")
	cmd.hide('lines', 'aln_seq2')
	cmd.hide('spheres', 'aln_seq2')
	cmd.hide('sticks', 'aln_seq2')
	

	cmd.hide('lines', 'star3d_aln')
	cmd.hide('spheres', 'star3d_aln')
	cmd.hide('sticks', 'star3d_aln')
	cmd.hide('ribbon', 'star3d_aln')

	cmd.color('red', 'aln_seq1')
	cmd.color('red', 'aln_seq2')
	
	for i in aln_index1:
		cmd.select('res', 'resi %s and chain A' % (i))
		cmd.color('yellow', 'res')

	for i in aln_index2:
		cmd.select('res', 'resi %s and chain B' % (i))
		cmd.color('blue', 'res')

	cmd.deselect()
	cmd.delete('res')

	cmd.bg_color('grey')	
cmd.extend("ShowAln", ShowAln)
