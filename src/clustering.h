#include <sys/stat.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>

#include "molecule.h"
#include "parameters.h"
#include "util.h"

#ifndef _HYDROCLUSTER_H
#define _HYDROCLUSTER_H

#define EPSILON 2

void WaterOccupancy(char *, ligand *, grid *, fgrid *, water *);

void Cluster_QT(ligand *, grid *, water *, hcluster *);
void QTClustering(grid *, hcluster *, water *);
void OutputWaterIndex(hcluster *, water *);

void Cluster_DB(fgrid *fgrd, ligand *lig, grid *grd, water *wtr, hcluster *cluster_head);
void dbScan(fgrid *fgrd, grid *grd, water *wtr, hcluster *cluster_head);
int  regionQuery(int , int,  int , int , float **, int *, int);
void expandCluster(float , int , int,  int , int , float **, int *, int *, int *, int, int *);
int  check(float *, int);

void Grid_cluster(grid *grd, water *wtr, hcluster *cluster_head);

void BulkCluster(fgrid *fgrd, grid *grd, water *wtr, hcluster *cluster_head);

#endif
