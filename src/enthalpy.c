#include "enthalpy.h"


//=================================================================================================
// CalculateEnthalpy
//=================================================================================================
void CalculateFreeEnergy(water *wtr, hcluster *cluster_0, energy *engy)
{

	hcluster	*cluster;
	float		deltaE, deltaS;
	float		enthalpy, TdeltaS;
	int		    *clwtrn;
	int         i, cl, frm, watindex;
	int         clusternum = cluster_0->numclusters;
	par.numclusters = cluster_0->numclusters;

	engy->enthalpy    = (float*)calloc(clusternum, sizeof(float));
	engy->occupancy   = (float*)calloc(clusternum, sizeof(float));
	engy->free_energy = (float*)calloc(clusternum, sizeof(float));
	clwtrn            = (int*)calloc(clusternum, sizeof(int));

	for(i = 0; i < clusternum; i++)
	{
		engy->enthalpy[i] = 0.0;
		clwtrn[i] = 0;
	}

	for(cluster = cluster_0, cl=0; cl < cluster_0->numclusters; cluster = cluster->next, cl++)
	{
		//printf("cluster %d \n", cl);
		for(frm = 0; frm < wtr->numfr; frm++)
		{
			if(cluster->wtrindex[frm] >= 0)
			{
				//printf("frame %d -- wtrindex %d\n", frm, cluster->wtrindex[frm] );
				watindex = cluster->wtrindex[frm];
				engy->enthalpy[cl] += wtr->dH[frm][watindex];
				clwtrn[cl]++;
				//printf("frame(%d): enthalpy=%f\n",frm, wtr->dH[frm][watindex]);
			}
		}
		engy->occupancy[cl] = clwtrn[cl] * 1.0 / wtr->numfr;
		//if(par.grid_energy == 0) printf("%d/%d = %f     ", clwtrn[cl], wtr->numfr, engy->occupancy[cl]);
		//printf("%d/%d = %f     \n", clwtrn[cl], wtr->numfr, engy->occupancy[cl]);
		engy->enthalpy[cl] /= clwtrn[cl];         //4.184: joules -> cal

		enthalpy = engy->enthalpy[cl];	//the cutoff is from Freisner's paper
		deltaE = enthalpy - par.E_bulk;

		TdeltaS = engy->entropy[cl];
		engy->enthalpy[cl] = deltaE;
		//deltaS = -TS;					//S_bulk - S_prot
		engy->free_energy[cl] = deltaE - TdeltaS;	//
		//if(par.grid_energy == 0) printf("cluster=%d :  dH (%f)  - TdS (%f) = dG (%f)\n", cl, deltaE, TdeltaS, engy->free_energy[cl]);
		//printf("cluster=%d :  dH (%f)  - TdS (%f) = dG (%f)\n", cl, deltaE, TdeltaS, engy->free_energy[cl]);
	}
}

//=================================================================================================
// Read SAMGP files for grid points
//=================================================================================================
void OutputClusterInfo(grid *grd, energy *engy)
{
	FILE		*fo, *foH, *foS, *foG;
	int		clustern;
	int		i, j;

	if(par.grid_energy == 1) // grid energy
	{

		foH = fopen("gridEgy_dH.dx", "w");
		foS = fopen("gridEgy_dTS.dx", "w");
		foG = fopen("gridEgy_dG.dx", "w");

		fprintf(foH,
				"object 1 class gridpositions counts %d %d %d\n"
				"origin %f %f %f\n"
				"delta %.2f 0.00 0.00\ndelta 0.00 %.2f 0.00\ndelta 0.00 0.00 %.2f\n"
				"object 2 class gridconnections counts %d %d %d\n"
				"object 3 class array type double rank 0 items %d data follows\n",
				grd->NumGP_r[0], grd->NumGP_r[1], grd->NumGP_r[2],
				grd->ZeroPoint[0], grd->ZeroPoint[1], grd->ZeroPoint[2],
				par.griddelta, par.griddelta, par.griddelta,
				grd->NumGP_r[0], grd->NumGP_r[1], grd->NumGP_r[2],
				grd->NumGP_r[0]*grd->NumGP_r[1]*grd->NumGP_r[2]);
		fprintf(foS,
				"object 1 class gridpositions counts %d %d %d\n"
				"origin %f %f %f\n"
				"delta %.2f 0.00 0.00\ndelta 0.00 %.2f 0.00\ndelta 0.00 0.00 %.2f\n"
				"object 2 class gridconnections counts %d %d %d\n"
				"object 3 class array type double rank 0 items %d data follows\n",
				grd->NumGP_r[0], grd->NumGP_r[1], grd->NumGP_r[2],
				grd->ZeroPoint[0], grd->ZeroPoint[1], grd->ZeroPoint[2],
				par.griddelta, par.griddelta, par.griddelta,
				grd->NumGP_r[0], grd->NumGP_r[1], grd->NumGP_r[2],
				grd->NumGP_r[0]*grd->NumGP_r[1]*grd->NumGP_r[2]);
		fprintf(foG,
				"object 1 class gridpositions counts %d %d %d\n"
				"origin %f %f %f\n"
				"delta %.2f 0.00 0.00\ndelta 0.00 %.2f 0.00\ndelta 0.00 0.00 %.2f\n"
				"object 2 class gridconnections counts %d %d %d\n"
				"object 3 class array type double rank 0 items %d data follows\n",
				grd->NumGP_r[0], grd->NumGP_r[1], grd->NumGP_r[2],
				grd->ZeroPoint[0], grd->ZeroPoint[1], grd->ZeroPoint[2],
				par.griddelta, par.griddelta, par.griddelta,
				grd->NumGP_r[0], grd->NumGP_r[1], grd->NumGP_r[2],
				grd->NumGP_r[0]*grd->NumGP_r[1]*grd->NumGP_r[2]);

		int counter;
		for (i = 0; i < grd->numGP3; i++) {
			counter = 0;
			//printf("i=%d\n", i);
			for(clustern = 0; clustern < par.numclusters; clustern++)
			{
				if (engy->grid_num[clustern] == i)
				{
					counter ++;
					//printf("cluster=%d, cluter_grid=%d, grd->numGP3=%d\n", clustern, engy->grid_num[clustern], i);
					fprintf(foH, "%10.5f",engy->enthalpy[clustern]);
					fprintf(foS, "%10.5f",engy->entropy[clustern]);
					fprintf(foG, "%10.5f",engy->free_energy[clustern]);
					break;
				}
			}
			if (counter == 0) {
				fprintf(foH, "%10.5f",0.0);
				fprintf(foS, "%10.5f",0.0);
				fprintf(foG, "%10.5f",0.0);
			}
			if ( (i+1)%3 == 0) {
				fprintf(foH, "\n");
				fprintf(foS, "\n");
				fprintf(foG, "\n");
			}
		}
		fclose(foH);
		fclose(foS);
		fclose(foG);
	

		foH = fopen("gridVis_dH.dx", "w");
		foS = fopen("gridVis_dTS.dx", "w");
		foG = fopen("gridVis_dG.dx", "w");

		fprintf(foH,
				"object 1 class gridpositions counts %d %d %d\n"
				"origin %f %f %f\n"
				"delta %.2f 0.00 0.00\ndelta 0.00 %.2f 0.00\ndelta 0.00 0.00 %.2f\n"
				"object 2 class gridconnections counts %d %d %d\n"
				"object 3 class array type double rank 0 items %d data follows\n",
				grd->NumGP_r[0], grd->NumGP_r[1], grd->NumGP_r[2],
				grd->ZeroPoint[0], grd->ZeroPoint[1], grd->ZeroPoint[2],
				par.griddelta, par.griddelta, par.griddelta,
				grd->NumGP_r[0], grd->NumGP_r[1], grd->NumGP_r[2],
				grd->NumGP_r[0]*grd->NumGP_r[1]*grd->NumGP_r[2]);
		fprintf(foS,
				"object 1 class gridpositions counts %d %d %d\n"
				"origin %f %f %f\n"
				"delta %.2f 0.00 0.00\ndelta 0.00 %.2f 0.00\ndelta 0.00 0.00 %.2f\n"
				"object 2 class gridconnections counts %d %d %d\n"
				"object 3 class array type double rank 0 items %d data follows\n",
				grd->NumGP_r[0], grd->NumGP_r[1], grd->NumGP_r[2],
				grd->ZeroPoint[0], grd->ZeroPoint[1], grd->ZeroPoint[2],
				par.griddelta, par.griddelta, par.griddelta,
				grd->NumGP_r[0], grd->NumGP_r[1], grd->NumGP_r[2],
				grd->NumGP_r[0]*grd->NumGP_r[1]*grd->NumGP_r[2]);
		fprintf(foG,
				"object 1 class gridpositions counts %d %d %d\n"
				"origin %f %f %f\n"
				"delta %.2f 0.00 0.00\ndelta 0.00 %.2f 0.00\ndelta 0.00 0.00 %.2f\n"
				"object 2 class gridconnections counts %d %d %d\n"
				"object 3 class array type double rank 0 items %d data follows\n",
				grd->NumGP_r[0], grd->NumGP_r[1], grd->NumGP_r[2],
				grd->ZeroPoint[0], grd->ZeroPoint[1], grd->ZeroPoint[2],
				par.griddelta, par.griddelta, par.griddelta,
				grd->NumGP_r[0], grd->NumGP_r[1], grd->NumGP_r[2],
				grd->NumGP_r[0]*grd->NumGP_r[1]*grd->NumGP_r[2]);

		float maxH = 15.0;
		float maxS = 15.0;
		for (i = 0; i < grd->numGP3; i++) {
			counter = 0;
			//printf("i=%d\n", i);
			for(clustern = 0; clustern < par.numclusters; clustern++)
			{
				if (engy->grid_num[clustern] == i)
				{
					counter ++;
					//printf("cluster=%d, cluter_grid=%d, grd->numGP3=%d\n", clustern, engy->grid_num[clustern], i);
					fprintf(foH, "%10.5f",engy->enthalpy[clustern]);
					fprintf(foS, "%10.5f",engy->entropy[clustern]);
					fprintf(foG, "%10.5f",engy->free_energy[clustern]);
					break;
				}
			}
			if (counter == 0) {
				fprintf(foH, "%10.5f",maxH);
				fprintf(foS, "%10.5f",maxS);
				fprintf(foG, "%10.5f",maxH+maxS);
			}
			if ( (i+1)%3 == 0) {
				fprintf(foH, "\n");
				fprintf(foS, "\n");
				fprintf(foG, "\n");
			}
		}
		fclose(foH);
		fclose(foS);
		fclose(foG);


	}
	else{

		// hydration site analysis output

		fo = fopen("cluster.egy", "w");

		fprintf(fo, "@cluster# entropy(-TS_e) enthalpy total(E-TS) center_x center_y center_z radius occupancy\n");
		for(clustern = 0; clustern < par.numclusters; clustern++)
		{
			fprintf(fo, "%4d %10.3f %10.3f %10.3f %7.4f %7.4f %7.4f %7.2f   %5.5f   %5.5f\n", clustern, -engy->entropy[clustern], engy->enthalpy[clustern], engy->free_energy[clustern], engy->coor[clustern][0], engy->coor[clustern][1], engy->coor[clustern][2], engy->radius[clustern], engy->occupancy[clustern], engy->overall_occupancy[clustern]);

			//printf("%d %f\n", clustern, engy->entropy[clustern]);
		}
		fclose(fo);

		fo = fopen("HydrationSites.pdb", "w");
		fprintf(fo, "@cluster# entropy(-TS_e) enthalpy total(E-TS) center_x center_y center_z radius occupancy\n");
		for(clustern = 0; clustern < par.numclusters; clustern++)
		{
			fprintf(fo, "ATOM    %3d  O   WAT A %3d    %8.3f%8.3f%8.3f%6.2f%6.2f           O\n", clustern+1, clustern+1, engy->coor[clustern][0], engy->coor[clustern][1], engy->coor[clustern][2], -engy->entropy[clustern], engy->free_energy[clustern]);
		}
		fclose(fo);

		fo = fopen("HydrationSites.mol2", "w");
		fprintf(fo, "@<TRIPOS>MOLECULE\n****\n");
		fprintf(fo, "% 5d % 5d     0     0     0\n", par.numclusters, 0);
		fprintf(fo, "SMALL\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n");
		for(clustern = 0; clustern < par.numclusters; clustern++)
		{
			fprintf(fo, "% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s %3d WAT       % 8.4f %5.2f\n",
					clustern+1, "WAT", engy->coor[clustern][0], engy->coor[clustern][1], engy->coor[clustern][2],
					"O", clustern+1, engy->free_energy[clustern], -engy->entropy[clustern]);
		}
		fclose(fo);

		fo = fopen("HS_dG.pdb", "w");
		fprintf(fo, "@cluster# entropy(-TS_e) enthalpy total(E-TS) center_x center_y center_z dH dG\n");
		for(clustern = 0; clustern < par.numclusters; clustern++)
		{
			fprintf(fo, "ATOM    %3d  O   WAT A %3d    %8.3f%8.3f%8.3f%6.2f%6.2f           O\n", clustern+1, clustern+1, engy->coor[clustern][0], engy->coor[clustern][1], engy->coor[clustern][2], engy->free_energy[clustern], engy->free_energy[clustern]);
		}
		fclose(fo);

		fo = fopen("HS_dH.pdb", "w");
		fprintf(fo, "@cluster# entropy(-TS_e) enthalpy total(E-TS) center_x center_y center_z dTS dH\n");
		for(clustern = 0; clustern < par.numclusters; clustern++)
		{
			fprintf(fo, "ATOM    %3d  O   WAT A %3d    %8.3f%8.3f%8.3f%6.2f%6.2f           O\n", clustern+1, clustern+1, engy->coor[clustern][0], engy->coor[clustern][1], engy->coor[clustern][2], -engy->entropy[clustern], engy->enthalpy[clustern]);
		}
		fclose(fo);

		fo = fopen("HS_dTS.pdb", "w");
		fprintf(fo, "@cluster# entropy(-TS_e) enthalpy total(E-TS) center_x center_y center_z dH dTS\n");
		for(clustern = 0; clustern < par.numclusters; clustern++)
		{
			fprintf(fo, "ATOM    %3d  O   WAT A %3d    %8.3f%8.3f%8.3f%6.2f%6.2f           O\n", clustern+1, clustern+1, engy->coor[clustern][0], engy->coor[clustern][1], engy->coor[clustern][2], engy->enthalpy[clustern], -engy->entropy[clustern]);
		}
		fclose(fo);


	}
}



