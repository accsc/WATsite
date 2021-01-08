#include <stdlib.h>
#include <stdio.h>

#ifndef _PARAMETERS_H
#define _PARAMETERS_H

#define kString 254
#define lString 508

#define Null 0.00001
#define Pi 				3.141592654
#define sqrtPi			1.772453851
#define sqrt_2			1.414213562
#define E 				2.718281828
// gas constant in cal/(K*mol)
#define R 				1.9858775		
// room temperature in kelvin
#define T0				298.15	
// degrees to radians		
#define	DtoR			0.0174532925199	
#define RtoD			57.2957795131

//enthalpy for bulk solvent water in kcal/mol (-8.13 from experimental data, -18.5 from Freisner paper)
#define BulkSolE		-8.13
//entropy for bulk solvent in cal/mol-k (16.42 from experimental data, 5.03 from Freisner paper)
#define BulkSolS		16.42

/*
#define BulkE_spce     0.00 
#define BulkE_tip3p    0.00 
#define BulkE_opc      0.00 
#define BulkE_tip4p    0.00 
#define BulkE_tip4pew  0.00 
#define BulkE_tip5p    0.00 
#define BulkE_amoeba   0.00
#define BulkS          0.00 
*/


#define BulkE_spce      -21.15
#define BulkE_tip3p     -17.95
#define BulkE_opc       -23.14
#define BulkE_tip4p     -18.70
#define BulkE_tip4pew   -20.95
#define BulkE_amoeba    -18.46
//#define BulkS           -3.78 from MD for HS
#define BulkS           -4.598
#define BulkE_tip5p      0.00 


struct parameters
{
    char		ParameterFile[lString];
    char		outfolder[lString];

    char		LigandFile[lString];
    char		ProtFile[lString];
    char		WaterEnthalpyFolder[lString];

	int		numfr;
	char		water_model[lString];
	int		site3_water, site4_water, site5_water;
	float	E_bulk;

	float		size_water;
	float		griddelta;
	float		FarestDist;
	float		griddensity;

	int			covdimension;

	int			cluster_method;
	int			grid_energy;

	int			dbscan_START;
	int			dbscan_END;

	int			number_clusters;
	float		cluster_mean_dist;
	float		maxdist;
	float		QT_occupancy_cutoff;

	int			occupancy_cutoff;
	int			numbins;
	//float		NeatS;

	//float		Ecutoff, Scutoff, Ereward, TSreward;
	//float		E_bulk;

	int			numwaters;
	int			numclusters;
	int			numframes;

	//MM 1/11/2018 Clustered trajectory parameters
	int			clusteredTrajFlag;				//set to 1 if the trajectory is clustered
	int			frameCutoffFlag;				//set to 1 if the first frames are cut
    char        TrajEgyAssocFile[lString];
	int			numclust;
	int			totalnumfr;
	int			frameCutoff;

} par;

#endif
