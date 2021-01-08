#!/usr/bin/python
#title           :truncate_prot.py
#description     :This will truncate protein by some distance around the binding site and cap terminals with ACE/NME.
#                 Any residue within such distance will be kept intact.
#author          :Ying Yang
#date            :11-13-2018
#usage           :python truncate_prot.py --pdb_amber prot_dry.pdb --bs MyBindingSite.pdb -t 15.0
#python_version  :2.7  
#==============================================================================

import os,sys
import numpy as np
from argparse import ArgumentParser
from collections import defaultdict
import pymol

pymol.pymol_argv = ['pymol','-qc'] #+ sys.argv[1:]
pymol.finish_launching()
cmd = pymol.cmd


parser = ArgumentParser(description='This will shift the geometry center of the system to the center of box.')
parser.add_argument('--pdb_propka',     dest='pdb_propka', metavar='<pdb FILE>', help='Input pdb file with propka protonation done. ')
parser.add_argument('--pdb_amber',      dest='pdb_amber',  metavar='<pdb FILE>', help='Input pdb file with amber f.f. compatible residue names. ')
parser.add_argument('--bs',             dest='bs', metavar='<pdb FILE>', help='Input pdb file to define binding site for truncation.')
parser.add_argument('-t', '--truncate', dest='truncate', metavar='<distance cutoff>',type=float, help='truncate protein by <input> distance cutoff via pymol')



def main():
	opt = parser.parse_args()
	if opt.pdb_propka:
		rename_chargeRes(opt.pdb_propka)
	if opt.truncate:
		truncate(opt.truncate)

def rename_chargeRes(filename):
	residues = ['ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'CYS', 'CYX', 'CY1', 'GLH', 'GLN', 'GLU', 'GLY', 'HIS', 'HIE', 'HID', 'HIP', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HEM', 'HEO']
	ions = ["NA", "Na+", "CL", "Cl-", "CA", "ZN"]
	wat = ["HOH", "WAT", "SOL"]

	prot_lines = []
	Fo_lig = open("tmp_lig.pdb", "w")
	Fo_ion = open("tmp_ion.pdb", "w")
	Fo_wat = open("tmp_wat.pdb", "w")
	with open(filename,"r") as Fi:
		for line in Fi:
			if line.find("ATOM") > -1 or line.find("HETATM") > -1:
				if line[17:20] in residues:
					prot_lines.append(line)
				elif line[17:20].strip() in ions:
					Fo_ion.write(line)
				elif line[17:20].strip() in wat:
					if line[12:16].strip()[0] != "H":
						Fo_wat.write(line)
				else:
					new_line = line[:17] + "SUB" + line[20:]
					Fo_lig.write(new_line)
	Fo_lig.close()
	Fo_ion.close()
	Fo_wat.close()

	out = []
	chain_list = []
	for line in prot_lines:
		if not line[21:22].upper() in chain_list:
			chain_list.append(line[21:22].upper())
	print "Protein pdb contains chain(s)", chain_list

	for chain in chain_list:
		his_dict = defaultdict(list)
		asp_dict = defaultdict(list)
		glu_dict = defaultdict(list)
		lys_dict = defaultdict(list)
		for line in prot_lines:
			if (line[0:4] == 'ATOM' or line[0:6] == 'HETATM') and line[21:22].upper() == chain :
				res_id = str( int(line[22:26]) )
				atom_name = line[12:16].strip()
				if line[17:20] == 'HIS':
					his_dict[res_id].append(atom_name)
				elif line[17:20] == 'ASP':
					asp_dict[res_id].append(atom_name)
				elif line[17:20] == 'GLU':
					glu_dict[res_id].append(atom_name)
				elif line[17:20] == 'LYS':
					lys_dict[res_id].append(atom_name)

		protna = {}
		for key in his_dict:
			if "HD1" in his_dict[key] and "HE2" in his_dict[key]:
				protna[key] = "HIP"
			elif "HD1" in his_dict[key] and not "HE2" in his_dict[key]:
				protna[key] = "HID"
			elif "HE2" in his_dict[key] and not "HD1" in his_dict[key]:
				protna[key] = "HIE"
			else:
				protna[key] = "HIM"
		for key in asp_dict:
			if "HD2" in asp_dict[key]:
				protna[key] = "ASH"
			else:
				protna[key] = "ASP"
		for key in glu_dict:
			if "HE2" in glu_dict[key]:
				protna[key] = "GLH"
			else:
				protna[key] = "GLU"
		for key in lys_dict:
			if "HZ3" in lys_dict[key]:
				protna[key] = "LYS"
			else:
				protna[key] = "LYN"
		print protna

		for line in prot_lines:
			if (line[0:4] == 'ATOM' or line[0:6] == 'HETATM') and line[21:22].upper() == chain :
				res_id = str(int(line[22:26]))
				atom_name = line[12:16].strip()
				alternative_loc = str(line[16])
				if atom_name[0] != "H" and alternative_loc != "B":
					if res_id in protna:
						out.append( line[0:17] + protna[res_id] + " " + line[21:] )
					else:
						out.append( line )
		out.append("TER\n")
	print 'HIS/ ASP/ GLU/ LYS residues are renamed'

	Fo_prot = open("tmp_prot.pdb", "w")
	for line in out:
		Fo_prot.write(line)
	Fo_prot.close()

	os.system("cat tmp_prot.pdb tmp_wat.pdb > prot.pdb")
	if os.stat("tmp_lig.pdb").st_size == 0:
		os.system("rm tmp_lig.pdb")
	else:
		os.system("mv tmp_lig.pdb bs_ligand.pdb")
	os.system("rm tmp_ion.pdb tmp_wat.pdb tmp_prot.pdb")

def truncate(dis):
	
	lines = []
	chain = 0
	chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	with open("protH.pdb") as protH:
		for line in protH:
			if line.split()[0] == "TER":
				chain += 1
			elif line.split()[0] == "ATOM" and line.split()[3] != "WAT":
				line = line[:21] + chains[chain] + line[22:]
			lines.append(line)

	with open("protH.pdb", "w") as out:
		for line in lines:
			out.write(line)

	cmd.load('protH.pdb')
	cmd.remove('(hydro)')
	cmd.load('MyBindingSite.pdb')
	cmd.select('br. MyBindingSite around %f and protH' % dis)
	cmd.save('prot_truncated.pdb', 'sele')

	inpdb = 'prot_truncated.pdb'
	fullpdb = 'protH.pdb'
	outpdb = "prot_prep.pdb"

	info  = {}
	lines = defaultdict(list)
	with open(fullpdb, "r") as f:
		for line in f:
			if line.find("ATOM") > -1 or line.find("HETATM") > -1:
				resNum, ele, x, y, z = int(line[22:26]), line[12:16].strip(), float(line[30:38]), float(line[38:46]), float(line[46:54])
				if ele in ['C','CA','O','N']:
					info['%d-%s' %(resNum, ele)] = np.array([x, y, z])
				lines[resNum].append(line)

	fi = open(inpdb, "r")
	fo = open(outpdb, "w")
	resNumList = []
	for line in fi:
		if (line.find("ATOM") > -1 or line.find("HETATM") > -1) and not (line.find("WAT") > -1):
			resNum, ele = int(line[22:26]), line[12:16].strip()
			if resNum not in resNumList: # N atom of a new residue
				resNumList.append(resNum)
				if len(resNumList) == 1: # first residue
					find_CA = "%d-CA" % (resNum-1)
					find_C  = "%d-C"  % (resNum-1)
					find_O  = "%d-O"  % (resNum-1)
					if find_CA in info:
						#print(resNum, resNum-1, lines[resNum - 1][0])
						#print(lines[ resNum - 1 ][0][21:22])
						chain_id = lines[ resNum - 1 ][0][21:22]

						fo.write("ATOM  %5d  CH3 ACE %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-1], info[find_CA][0], info[find_CA][1], info[find_CA][2]) )
						fo.write("ATOM  %5d  C   ACE %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-1], info[find_C][0], info[find_C][1], info[find_C][2]) )
						fo.write("ATOM  %5d  O   ACE %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-1], info[find_O][0], info[find_O][1], info[find_O][2]) )
					else:
						### MAY CAUSE PROBLEM
						pass
						# don't write at the beginning of original terminal # YY 3-27
						#fo.write("ATOM  %5d  C   ACE  %4d    %8.3f%8.3f%8.3f\n" % (0, resNumList[-1], 0, 0, 0) )
					fo.write(line)
				else:
					if resNumList[-1] == resNumList[-2] + 2: # only one residue deleted, add it back --> we don't break 
						res_add = resNumList[-2] + 1
						for line in lines[res_add]:
							fo.write(line)
						fo.write(line)
					elif resNumList[-1] == resNumList[-2] + 1: # no break
						fo.write(line)
					else: # this resnum is larger than prev +1 
						find_N  = "%d-N"  % (resNumList[-2] + 1)
						find_CA = "%d-CA" % (resNumList[-2] + 1)
						chain_id = lines[ resNumList[-2] + 1 ][0][21:22]
						fo.write("ATOM  %5d    N NME %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-2], info[find_N ][0], info[find_N ][1], info[find_N ][2]) )
						fo.write("ATOM  %5d  CH3 NME %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-2], info[find_CA][0], info[find_CA][1], info[find_CA][2]) ) 

						fo.write("TER\n")

						find_CA = "%d-CA" % (resNumList[-1]-1)
						find_C  = "%d-C"  % (resNumList[-1]-1)
						find_O  = "%d-O"  % (resNumList[-1]-1)
						chain_id = lines[ resNumList[-1]-1 ][0][21:22]
						fo.write("ATOM  %5d  CH3 ACE %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-1], info[find_CA][0], info[find_CA][1], info[find_CA][2]) )
						fo.write("ATOM  %5d  C   ACE %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-1], info[find_C ][0], info[find_C ][1], info[find_C ][2]) )
						fo.write("ATOM  %5d  O   ACE %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-1], info[find_O ][0], info[find_O ][1], info[find_O ][2]) )

						fo.write(line)
			else: # for other atoms in the same residue
				fo.write(line)
		elif line.find("END") > -1:
			resNum = resNumList[-1] # prev res
			find_N  = "%d-N"  % (resNum+1)
			find_CA = "%d-CA" % (resNum+1)
			if find_CA in info:
				chain_id = lines[ resNum+1 ][0][21:22]
				fo.write("ATOM  %5d    N NME %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-1], info[find_N ][0], info[find_N ][1], info[find_N ][2]) )
				fo.write("ATOM  %5d  CH3 NME %s%4d    %8.3f%8.3f%8.3f\n" % (0, chain_id, resNumList[-1], info[find_CA][0], info[find_CA][1], info[find_CA][2]) ) 
			else:
				pass
			fo.write("END\n")
		elif line.find("TER") > -1:
			pass
		else: # not ATOM/HETATM
			fo.write(line)
	fi.close()
	fo.close()

	fo = open('tmp', "w")
	ter = 0
	with open(outpdb, "r") as fi:
		for line in fi:
			#if (line.find("OXT") > -1 or line.find(" H2 ") > -1 or line.find(" H3 ") > -1 ) and not line.find("WAT") > -1:
			# OXT or H1-3 should not be added if tleap runs correctly. // YY 3-28
			# It was added because not chain was specified in the cap.pdb file
			"""
			if not line.find("WAT") > -1:
				pass
			else:
				if line.find("WAT") > -1 and ter==0:
					fo.write("TER\n")
					fo.write(line)
					ter = 1
				else:
					fo.write(line)
			"""
			if line.find("WAT") > -1 and ter==0:
				fo.write("TER\n")
				fo.write(line)
				ter = 1
			else:
				fo.write(line)

	fo.close()
	os.system("mv tmp %s" % outpdb)
	#os.system("sed -i /OXT/d %s" % outpdb)
	


def gen_files():

	with open('sys_prep', "w") as f:
		f.write("3.0000\n")
		f.write("amber14SB\n")
		f.write("SPC/E\n")

	with open('WATsite.out', "w") as f:
		f.write("ProductionStep      10000000\n")
		f.write("Interval            500\n")
		f.write("numFrame            20000\n")
		f.write("Protein             MyProtein.pdb\n")
		f.write("HydrationSiteMol2   WATsite_OUT/HydrationSites.mol2\n")
		f.write("HydrationSitePDB    WATsite_OUT/HydrationSites.pdb\n")
		f.write("GridMol2            WATsite_OUT/FilterGrid.mol2\n")
		f.write("GridDX              WATsite_OUT/grid_occupancy.dx\n")
		f.write("WaterTraj           WATsite_OUT/WATinside.mol2\n")
		f.write("EnergyFile          WATsite_OUT/cluster.egy\n")
		f.write("\n")
		f.write("\n")

	os.system("cp bs_ligand.pdb MyBindingSite.pdb")
	os.system("cp truncated_cap.pdb MyProtein.pdb")


if __name__ == "__main__":
    main()

