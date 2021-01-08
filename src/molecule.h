#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "parameters.h"

#ifndef _MOLECULE_H
#define _MOLECULE_H

typedef struct atom {
	int Num; // atom number start from 0
	char name[8];
	float x[3];
	int nearestgrdp;  //nearest grid point of this atom
} atom;

typedef struct ligand {
	atom *atm;
	int numheavyatms;
} ligand;

typedef struct grid {
	float xmin[3], xmax[3];
	float griddelta;     // griddelta: distance between two adjacent grid points
	float GridSize[3];             // grd->xmax[j] - grd->xmin[j];
	float ZeroPoint[3];            // grd->ZeroPoint[j] = grd->xmin[j];

	int   NumGP_r[3];
	int   numGP1, numGP2, numGP3;	// grd->numGP3 = grd->NumGP_r[0]*grd->NumGP_r[1]*grd->NumGP_r[2];

	float **GP_r;                  // store the coordinates of grid points
								   // grd->GP_r[num][0-2] x,y,z coordinates of the num_th grid point
	float **interaction;

	int   *flag;
	float *WaterOccupancy;      // Water Occupancy of each grid points
	float *Entropy;             // Water Entropy of each grid points
	float *Enthalpy;            // Water Enthalpy of each grid points
} grid;

typedef struct water//this structure only have the information of waters inside of binding pocket
{
	int numfr;	        //number of frames in snapshots
	int numatms;	    //number of total atoms for one frame in NMR-PDB
	int numWATs;//number of total waters for one frame (equals to number of water oxygen)

	atom **O;		    //oxygen coordinates for waters inside binding pocket
	atom **H1, **H2, **H3, **H4;//hydrogen coordinates for waters inside binding pocket (H3, H4: dummy atom for TIP4P, TIP5P waters)

	int **atn;//atom numbers of water inside binding pocket (corresponding to WATinside: 0, 1...)
	float **dH;	//enthalpy (deltaH) of water inside binding pocket (corresponding to WATinside: 0, 1...)
	int **resn;	//sol number (resn) of water inside binding pocket (for the purpose of debugging)
	int *WATinside;	//number of waters inside of binding pocket for each frame

	int largest_wat;//largest number of waters inside binding pocket among all the frames
	int **nearestgrdp;	//nearest grid point for each water in each frame
	int totalresn;
} water;

typedef struct filtered_grid {
	float **GP_r;
	float griddelta;
	int gridnum;
	int *flag;
	float *WaterOccupancy;
	int *ndx_in_grid;
} fgrid;

typedef struct hydro_cluster {
	int		num;
	int		numclusters;	//number of clusters in total
	float	**hydropoint;	//coor of grid point inside the cluster
	float	center[3];
	float	*hydroscore;	//occupancy for each point
	int		numpoints;
	float	average_occupancy;			//average occupancy of the cluster
	float	center_occupancy;	//occupancy of the center point
	float	max_occupancy;	//maximum occupancy among the cluster points
	float	overall_occupancy;
	int		*grdpn;			//grid point number
	int		grid_num;			//grid point number
	int		*wtrindex;//water number inside of clusters for each frame (number corresponding to order in water structure)
	int		wtrnum;			//number of waters inside of this cluster
	float		radius;
	float	**wtrEulerAgl;	//Euler angles

	float	cov[6][6];
	float	eigenvt[6][6];
	float	eigenvl[6];

	float	**PCAcoor;		//coordinate after PCA projection
	float	interval[6];	//bin size on each principle component dimension
	float	min[6], max[6];	//minimum and maximum in each PC dimension
	float	integral[6];
	float	probability[6][70];		//probability density of each bin
	struct	hydro_cluster *next;
} hcluster;

typedef struct energy {
//	int		    **atomlist;
//	int		    **resnlist;
	int		*grid_num;
	float *enthalpy;
	float *entropy;
	float *free_energy;
	float **coor;
	float *radius;
	float *occupancy;
	float *overall_occupancy;
	float *average_occupancy;
} energy;

#endif
