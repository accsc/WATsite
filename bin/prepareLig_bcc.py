#!/usr/bin/python
#title           :prepareLig_bcc.py
#description     :This will parameterize input ligand with GAFF using antechamber from AmberTools.
#author          :Ying Yang
#date            :11-13-2018
#usage           :python prepareLig_bcc.py -l MyLigand.mol2 -c 0
#python_version  :2.7  
#==============================================================================

import os, sys, re, getopt, glob
from argparse import ArgumentParser
import subprocess

parser = ArgumentParser(description='This will parameterize input ligand with GAFF using antechamber from AmberTools.')
parser.add_argument('-l','--lig', dest='lig', default='MyLigand.mol2', help='Input ligand mol2 file')
parser.add_argument('-c','--charge', dest='lig_charge',  help='Net charge of the input ligand')

def main():
    opt = parser.parse_args()

    if opt.lig:
        lig_file = opt.lig
        lig_charge = int(opt.lig_charge)
        out_dir = 'ligprep'

        print "Input ligand file: ", lig_file
        print "Input ligand chanrge: ", lig_charge
        print "Output will be generated in: ", out_dir
	
        try:
            os.mkdir(out_dir)
        except OSError as e:
            if e.errno == 17: # errno.EEXIST
                print("Warning: (%s) exists" % out_dir)
            else:
                raise

        fi = open(lig_file, "r")
        fo = open("ligprep/lig.mol2", "w")

        line = fi.readline()
        while line != '':
            if line.find("@<TRIPOS>MOLECULE") > -1:
                fo.write(line)
                line = fi.readline()
                fo.write(line)
                line = fi.readline()
                fo.write(line)
                numatoms = int(line.split()[0])
            elif line.find("@<TRIPOS>ATOM") > -1:
                fo.write(line)
                hcount = 1
                for i in range(numatoms):
                    line = fi.readline()
                    ss = line.split()
                    atom_name = ss[1]
                    newline = "% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  1 SUB       % 8.4f\n"%(int(ss[0]), atom_name, float(ss[2]), float(ss[3]), float(ss[4]), ss[5], float(ss[-1]))
                    fo.write(newline)
            else:
                fo.write(line)
            line = fi.readline()
        fi.close()
        fo.close()

        os.chdir(out_dir)
        genTleap(lig_file, lig_charge)
        os.chdir("../")

        command = 'cat prot.pdb ligprep/lig.pdb > complex.pdb'
        subprocess.call(command, shell=True)
        command = 'sed -i "/END/d" complex.pdb'
        subprocess.call(command, shell=True)


        fo = open('tmp.pdb', 'w')
        fi = open('complex.pdb', 'r')
        lines = fi.readlines()
        fi.close()
        for i in range( len(lines) ):
            if lines[i].find("OXT") > -1:
                fo.write(lines[i])
                fo.write("TER\n")
            elif lines[i].find("HOH") > -1 and not lines[i-1].find("HOH") > -1:
                fo.write("TER\n")
                fo.write(lines[i])
            elif lines[i].find("SUB") > -1 and not lines[i-1].find("SUB") > -1:
                fo.write("TER\n")
                fo.write(lines[i])
            else:
                fo.write(lines[i])
        fo.close()
        os.system("mv tmp.pdb complex.pdb")

# from rdkit import Chem
#def getNetCharge(lig_file):
#    m = Chem.MolFromMol2File(lig_file)
#    charge = Chem.GetFormalCharge(m)
#    return charge

def genTleap(lig_file, lig_charge):

    if os.path.isfile("lig_bcc_gaff.mol2"):
        print "\n\n"+"#"+"="*25+"#"
        print "Found output from previous antechamber run.\nDelete the previous project or create a new project if you want to perform antechamber again.\n"
    else:
        cm = '-c bcc'
        command = 'antechamber -i lig.mol2 -fi mol2 -o lig_bcc_gaff.mol2 -fo mol2 %s -nc %i -s 2 -at gaff -pf y -j 4' % (cm, lig_charge)
        print "Executing Antechamber..."
        print command
        subprocess.call(command, shell=True)

        command = 'antechamber -i lig.mol2 -fi mol2 -o lig.pdb -fo pdb -pf y '
        subprocess.call(command, shell=True)

    if not os.path.isfile("lig_bcc_gaff.mol2"):
        print "\nError: antechamber aborting ! "
        hint1 = "HINT1: is path of amber_home correctly set?\n"
        hint2 = "HINT2: is ligand charge correctly set?\n"
        print hint1,hint2
        sys.exit(1)

    if not os.path.isfile("lig_AC.frcmod"):
        command = 'parmchk2 -i lig_bcc_gaff.mol2 -f mol2 -s 1 -o lig_AC.frcmod' 
        print "Executing parmchk...\n %s" % command
        subprocess.call(command, shell=True)

if __name__ == "__main__":
    main()



