#!/usr/bin/python
#title           :get_centroid_pdb.py
#description     :This will compute the centroid of given pdb (binding site).
#author          :Ying Yang
#date            :11-13-2018
#usage           :python get_centroid_pdb.py MyBindingSite.pdb
#python_version  :2.7  
#==============================================================================


import sys, os
import numpy as np

def main(filename):
    coords = get_coords(filename)
    centroid_coords = centroid(coords)
    print("%.3f,%.3f,%.3f" % (centroid_coords[0], centroid_coords[1], centroid_coords[2]))
    
def get_coords(filename):
    fd = open(filename, "r")
    pdb_array = [ [line[30:38],line[38:46],line[46:54]] for line in fd.readlines() if line.startswith('ATOM')]
    #pdb_array = [line.split() for line in fd.readlines() if line.startswith('ATOM')]
    fd.close()

    coords = np.array([ map(float, atom[0:3]) for atom in pdb_array])
    #coords = np.array([ map(float, atom[6:9]) for atom in pdb_array])
    return coords

def centroid(atom_coords):
    centroid_coords = np.average(atom_coords, axis=0)
    assert centroid_coords.ndim == 1
    return centroid_coords

if __name__=='__main__':
    main(sys.argv[1])

