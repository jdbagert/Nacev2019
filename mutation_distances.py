############################################################################
############################################################################
#
#	mutation_distances.py
#
#	Calculates the distance between alpha carbons in highly mutated residues
#	in the nucleosome structure (1kx5.pdb) using Chimera (v 1.11.2) and 
#	writes the data to an output file
#
#	John D. Bagert
#	7/23/2018
#
############################################################################
############################################################################

import chimera
from chimera import runCommand as rc
import numpy as np
from matplotlib import cm

### Functions
############################################################################
### Calculates the distance between two atoms in Chimera
def atom_distance(spec1, spec2):
    atom1 = chimera.specifier.evalSpec(spec1).atoms()[0]
    atom2 = chimera.specifier.evalSpec(spec2).atoms()[0]
    d = atom1.xformCoord().distance(atom2.xformCoord())
    return d
############################################################################

### Script operations
writeoutput = 1 													# write output distances to file
inputfile = '18_07_13_oncohistone_Analysis_all_perpatient.txt'		# data input filename
outputfile = 'mut_dist_2.5x_10tmb_novar_notails.txt'				# data output filename

### User-defined parameters
mutcutoffmultiplier = 2.5			# mutational frequency cutoff multiplier (of median value for each histone) to cutoff mutants analyzed
distcutoff = 20.					# distance over which values are not recorded
tmbcutoff = 10.						# tumor mutational burden cutoff
exclude_h2a_variants = True			# excludes H2A variants
h2avars = ['P0C5Y9','Q71UI9','P0C0S5','P16104','O75367']	# Uniprot identifiers of H2A variants to exclude from analysis
exclude_tails = True				# excludes unstructured domains (for proximity analysis)
globular_domains = {'H2A':[16,116],'H2B':[37,122],'H3':[44,131],'H4':[24,98]}	# residues of the histone globular domains


### initializes data storage variables
############################################################################
f_seq = open('histone_seq.txt') 	# opens file with histone uniprot IDs and sequences
protseq = {}						# stores histone sequences
for l in f_seq:
	[_uniprot, _seq] = l.strip().lower().split('\t')
	protseq[_uniprot] = _seq

protmutcount = {'H3':[],'H4':[],'H2A':[],'H2B':[]}	# holds mutation counts (aligned) for each histone
for _hist in protmutcount:
	protmutcount[_hist] = [0]*len(protseq[_hist.lower()])

### opens mutations data file and stores information
############################################################################
f_mut = open(inputfile)
fh_mut = f_mut.readline().lower().strip().split('\t')
for l in f_mut:
	# Stores line as temporary variables
	d=dict(zip(fh_mut,l.lower().strip().split('\t')))
	_uniprot = d['uniprot'].upper()
	_aaposadj = int(d['mutation position aligned'])
	_histoneshort = d['histone_short'].upper()
	_len = len(protseq[_histoneshort.lower()])
	_tmb = float(d['tumor_mutational_burden (nonsynonymous mutations per mb)'])
	
	if _tmb > tmbcutoff: continue
	if exclude_h2a_variants:
		if _uniprot in h2avars: continue
	
	if _histoneshort == 'H2B': _aaposadj +=3 # Corrects for the -3 AA nubmering shift in H2B mutants of 1kx5.pdb
	if _aaposadj > _len: continue	# Skips mutants longer than canonical histone seq
	
	# stores total mutation count data
	protmutcount[_histoneshort][_aaposadj-1] += 1
f_mut.close()

### initializes variables for working in chimera
############################################################################
cutoffs = {} # mutation rate cutoff for distance calculations
for _hist in protmutcount: 
	cutoffs[_hist]=np.median(protmutcount[_hist])*mutcutoffmultiplier
	if cutoffs[_hist] == 0: cutoffs[_hist] = 2.5
print cutoffs

# dictionaries to convert between chains and histones and vice versa
chain2hist = {
'A':'H3-1','B':'H4-1','C':'H2A-1','D':'H2B-1',
'E':'H3-2','F':'H4-2','G':'H2A-2','H':'H2B-2',
'I':'DNA-1', 'J':'DNA-2'
}
hist2chain= {
'H3':['A','E'], 'H4':['B','F'], 'H2A':['C','G'], 'H2B':['D','H']
}
# conversion dictionary for AA naming
aa2res = {	
'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE',
'G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU',
'M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG',
'S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR',}

# initializes variables to hold all distance information
distances = {'H3':{},'H4':{},'H2A':{},'H2B':{}} # holds the distance information
distheads = {'H3':[],'H4':[],'H2A':[],'H2B':[]} # holds the header information for distances


### Chimera commands
############################################################################

### opens PDB molecule and stores atom, bond, residue info
pdb = chimera.openModels.open('1kx5.pdb')
mol = pdb[0]
atomlist = mol.atoms
bondlist = mol.bonds
reslist  = mol.residues

### Cosmetics of the nucleosome structure
rc('col cornflower blue :.A-H')
rc('~disp solvent')
rc('~disp ions')
rc('nuc side ladder :.I-J')

print 'Sites to include in distance analysis'
print 'Histone\tResidue\tMutations'
### Measures distances for each mutation
for _hist in protmutcount:
	_chains = hist2chain[_hist.upper()] # chain of current histone
	for _aa in range(len(protmutcount[_hist])): # current amino acid
		_muts = protmutcount[_hist][_aa] # current number of mutations
		if _muts < cutoffs[_hist]: continue # skips residues with mutations fewer than cutoff
		_aa += 1
		if _hist == 'H2B': _aapdb = _aa - 3 # corrects for -3 numbering in H2B for 1kx5 structure
		else: _aapdb = _aa
		if _aapdb < 1: continue
		
		if exclude_tails:	# exlucdes unstructured domains in histones
			[_beginglob,_endglob] = globular_domains[_hist] 
			if _aa < _beginglob: continue
			if _aa > _endglob: continue
		
		print _hist,'\t',_aa,'\t',_muts
		distheads[_hist].append(_aa)
		distances[_hist][str(_aa)] = {'A':{},'B':{},'C':{},'D':{},'E':{},'F':{},'G':{},'H':{}}
		_ca = ':%d.%s@ca' %(_aapdb,_chains[0]) #current residue alpha carbon position

		# for each AA loops through all other AAs and calculates distances
		for _hist2 in protmutcount:
			_chains2 = hist2chain[_hist2.upper()]	# second chain
			for _aa2 in range(len(protmutcount[_hist2])):	# second chain amino acid
				_muts2 = protmutcount[_hist2][_aa2]	# second number of mutations
				if _muts2 < cutoffs[_hist2]: continue # skips residues with mutations fewer than cutoff
				_aa2 += 1
				if _hist2 == 'H2B': _aapdb2 = _aa2 - 3	# corrects for -3 numbering in H2B for 1kx5 structure
				else: _aapdb2 = _aa2
				if _aapdb2 < 1: continue
				for _ch in _chains2:
					_ca2 = ':%d.%s@ca' %(_aapdb2,_ch) #current residue alpha carbon pos
					_dist = atom_distance(_ca,_ca2)
					distances[_hist][str(_aa)][_ch][str(_aa2)] = _dist

### Writes distance data to an output file
if writeoutput:
	fout = open(outputfile,'w')
	for _hist in distheads:
		fout.write('\t')
		#fout.write('\t'.join([str(x) for x in distheads[_hist]]))
		fout.write('\t'.join([_hist+' '+str(x) for x in distheads[_hist]])) # includes histone in label
	fout.write('\n')
	for _hist in distheads:
		for _aa in distheads[_hist]:
			#fout.write('%s\t' %_aa)
			fout.write('%s %s\t' %(_hist,_aa)) # includes histone in label
			for _hist2 in distheads:
				for _aa2 in distheads[_hist2]:
					_ch = hist2chain[_hist2.upper()]
					_dist = min([distances[_hist][str(_aa)][_ch[0]][str(_aa2)],distances[_hist][str(_aa)][_ch[1]][str(_aa2)]])
					fout.write('%.1f\t' %_dist)
			fout.write('\n')
	fout.close()



print 'DONE!!!'