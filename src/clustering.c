/***************************************************************************
 *   This script
 *   1) calculates water occupancy by going through the MD snapshot and    *
 *   mapping occupancy onto grid using gaussian distribution.              *
 *   2) clusters grip points by either modified k-mean cluster algorithm
 or DBSCAN algorithm *
 ***************************************************************************/
#define _FILE_OFFSET_BITS 64
#include "clustering.h"

//=================================================================================================
// WaterOccupancy // read in NMRsnapshots, calculate wateroccupancy on grid: nmrpdb in vmd format
// modified for the water distribution: radius = 1A = 3*sigma (more steep distribution than 
// WaterOccupancy2)
//=================================================================================================
void WaterOccupancy(char *infile, ligand *lig, grid *grd, fgrid *fgrd, water *wtr) {
	FILE *fi, *fo;
	char line[lString], atmn[8], coor[kString];
	int numfr = 0, numatm = 0, numWAT = 0, cc, dd, countAtm, fstWATatm; //YY 4.21.2015
	float x, y, z, alpha, beta, gamma, tx, ty, tz, delta, dis, dis2, dx, dy, dz, xx, yy, zz, dd2;
	int gx, gy, gz, num, gnum, del, atn, resn;
	int i, j, k, k2, l2, m2;
	int watinside = 0;
	int G_num, closestGP;
	float size_water = par.size_water;
	float size_water2 = pow(size_water, 2);
	float sigma = size_water / 3.0;
	float sigma2 = pow(sigma, 2);
	float denom;
	int grid_counter, new_counter;

	delta = grd->griddelta;     // grid spacing
	dis = size_water;// cutoff radius around center of water molecule (3*sigma will account for 99.7% of the set)
	del = (int) (dis / delta) + 1;// number of grid points around center of water molecule up to cutoff
	fi = fopen(infile, "r");
	printf("start reading protein snapshots...\n");

	if (!fi) {
		printf("protein file %s is missing.\n", infile);
		exit(0);
	} 
	else {
		grd->WaterOccupancy = (float*) calloc(grd->numGP3, sizeof(float));
		for (i = 0; i < grd->numGP3; i++) {
			grd->WaterOccupancy[i] = 0;
		}

		// get number of atoms and waters inside of binding pocket
		rewind(fi);
		while (!feof(fi)) {
			fgets(line, lString, fi);
			if (strstr(line, "ATOM")) {
				numatm++;
			}
			if (strstr(line, "O   WAT") || strstr(line, "OW")) {
				numWAT++;
				GetChars(line, 30, 37, coor);
				x = atof(coor);
				GetChars(line, 38, 45, coor);
				y = atof(coor);
				GetChars(line, 46, 54, coor);
				z = atof(coor);

				if ((grd->xmin[0] < x && x < grd->xmax[0])
						&& (grd->xmin[1] < y && y < grd->xmax[1])
						&& (grd->xmin[2] < z && z < grd->xmax[2]))
					watinside++;
			}
			if (strstr(line, "END"))
				break;
		}
		//wtr->totalresn = resn;
		wtr->numatms = numatm;
		wtr->numWATs = numWAT;
		printf("number of atoms:%d number of waters: %d\n", numatm, numWAT);
		printf("number of waters inside pocket in the first frame: %d\n", watinside);
		// get number of frames
		wtr->numfr = par.numfr;

		FILE *popen(const char *command, const char *mode);
		int pclose(FILE *stream);
		FILE *cmd;
		char result[1024];

		char string[100];
		sprintf(string, "tail -n 100000 %s | grep MODEL", infile);
		//printf("%s\n", string);
		cmd = popen(string, "r");
		if (cmd == NULL) {
			perror("popen");
			exit(EXIT_FAILURE);
		}
		int total_frames;
		while (fgets(result, sizeof(result), cmd)) {
			//printf("%s", result);
			sscanf(result, "%*s%d", &total_frames);
		}
		pclose(cmd);
		printf("Total number of frames from grep: %d\n", total_frames);

		if (total_frames > par.numfr) {
			printf("Warning: Number of frames in trajectory (%d) is larger than user input (%d)...\n", total_frames, par.numfr);
			numfr = par.numfr;
		} else if (total_frames < par.numfr) {
			printf("  Error: Number of frames in trajectory (%d) is smaller than user input (%d)...\n", total_frames, par.numfr);
			printf("Cannot proceed... Exit...\n");
			exit(0);
		} else{
			printf("User input frame number: %d \n", par.numfr);
			numfr = par.numfr;
		}
		rewind(fi);
		printf("number of frames %d\n", wtr->numfr);
		wtr->O = (atom**) calloc(numfr, sizeof(atom*));
		wtr->H1 = (atom**) calloc(numfr, sizeof(atom*));
		wtr->H2 = (atom**) calloc(numfr, sizeof(atom*));
		wtr->H3 = (atom**) calloc(numfr, sizeof(atom*));
		wtr->H4 = (atom**) calloc(numfr, sizeof(atom*));
		wtr->atn = (int**) calloc(numfr, sizeof(int*));
		wtr->resn = (int**) calloc(numfr, sizeof(int*));
		wtr->WATinside = (int*) calloc(numfr, sizeof(int));

		for (i = 0; i < numfr; i++) {
			wtr->O[i] = (atom*) calloc(watinside * 2, sizeof(atom));
			wtr->H1[i] = (atom*) calloc(watinside * 2, sizeof(atom));
			wtr->H2[i] = (atom*) calloc(watinside * 2, sizeof(atom));
			wtr->H3[i] = (atom*) calloc(watinside * 2, sizeof(atom));
			wtr->H4[i] = (atom*) calloc(watinside * 2, sizeof(atom));
			wtr->atn[i] = (int*) calloc(watinside * 2, sizeof(int));
			wtr->resn[i] = (int*) calloc(watinside * 2, sizeof(int));
		}

		// go over every frame of MD simulation, assign water occupancy to each grid point and store water coordinates
		for (i = 0; i < numfr; i++) {
			dd = 0;
			cc = 0;
			fstWATatm = 0; //YY 4.21.2015
			countAtm = 0; //YY 4.21.2015
			printf("%d / %d \r", i, numfr);
			while (!feof(fi)) {
				fgets(line, lString, fi);

				if (strstr(line, "O   WAT") || strstr(line, "OW  WAT")
						|| strstr(line, "O   HOH") || strstr(line, "OW  HOH")
						|| strstr(line, "OW  SOL")) {
					sscanf(line, "%*s%d%s%*s%d", &atn, atmn, &resn);
					if (countAtm == 0) //YY 4.21.2015
						fstWATatm = atn; //YY 4.21.2015
					//printf ("First Water Atom number is %d\n", fstWATatm); //YY 4.21.2015
					//resn = 0;
					GetChars(line, 30, 37, coor);
					x = atof(coor);
					GetChars(line, 38, 45, coor);
					y = atof(coor);
					GetChars(line, 46, 54, coor);
					z = atof(coor);

					// calculate occupancy for waters within defined binding pocket
					alpha = modff((x - grd->ZeroPoint[0]) / delta + Null, &tx);
					beta = modff((y - grd->ZeroPoint[1]) / delta + Null, &ty);
					gamma = modff((z - grd->ZeroPoint[2]) / delta + Null, &tz);
					//printf("Water oxygen atom: %f %f %f  %f %f %f  %.2f %.2f %.2f  %.2f %.2f %.2f\n",
					//      x, y, z, alpha, beta, gamma, tx, ty, tz, gx, gy, gz);
					gx = (int) (tx + Null);
					gy = (int) (ty + Null);
					gz = (int) (tz + Null);

					//YY grd->NumGP_r[0] = number of grid points on x axis

					//YY water outside binding pocket
					if (gx < 0 || gx >= grd->NumGP_r[0] - 1
					 || gy < 0 || gy >= grd->NumGP_r[1] - 1
					 || gz < 0 || gz >= grd->NumGP_r[2] - 1) 
					{
						cc++;
						// skip the next two hydrogen
						fgets(line, lString, fi);
						fgets(line, lString, fi);
						if (par.site4_water)
							fgets(line, lString, fi); //4 site water model
						if (par.site5_water)
							fgets(line, lString, fi); //5 site water model
					}

					else  //YY water inside binding pocket
					{   //YY go through the nearby grid points of this inside water  for each dimension x,y,z
						for (k2 = gx - del; k2 <= gx + del; k2++) {
							for (l2 = gy - del; l2 <= gy + del; l2++) {
								for (m2 = gz - del; m2 <= gz + del; m2++) {
									if (k2 >= 0 && k2 < grd->NumGP_r[0]
									 && l2 >= 0 && l2 < grd->NumGP_r[1]
									 && m2 >= 0 && m2 < grd->NumGP_r[2])
									{
										num = k2 * grd->numGP2 + l2 * grd->numGP1 + m2;
										xx = grd->GP_r[num][0] - x;
										yy = grd->GP_r[num][1] - y;
										zz = grd->GP_r[num][2] - z;
										dd2 = xx * xx + yy * yy + zz * zz;
										if (dd2 <= size_water2) {
											grd->WaterOccupancy[num] += exp(-dd2 / (2 * sigma2));
										}
									}
								}
							}
						}
						// assign nearest grid point to the water
						if (alpha < Null) {
							alpha = Null;
						}
						if (beta < Null) {
							beta = Null;
						}
						if (gamma < Null) {
							gamma = Null;
						}
						closestGP = (int) (2 * alpha) + 2 * (int) (2 * beta)
								+ 4 * (int) (2 * gamma); // this is like read, write, execute
						//printf("%f*2 = %f; int = %d\n", alpha, alpha*2, (int)(2*alpha) );
						switch (closestGP) {
						case 0:
							G_num = gx * grd->numGP2 + gy * grd->numGP1 + gz;
							break;
						case 1:
							G_num = (gx + 1) * grd->numGP2 + gy * grd->numGP1 + gz;
							break;
						case 2:
							G_num = gx * grd->numGP2 + (gy + 1) * grd->numGP1 + gz;
							break;
						case 3:
							G_num = (gx + 1) * grd->numGP2 + (gy + 1) * grd->numGP1 + gz;
							break;
						case 4:
							G_num = gx * grd->numGP2 + gy * grd->numGP1 + (gz + 1);
							break;
						case 5:
							G_num = (gx + 1) * grd->numGP2 + gy * grd->numGP1 + (gz + 1);
							break;
						case 6:
							G_num = gx * grd->numGP2 + (gy + 1) * grd->numGP1 + (gz + 1);
							break;
						case 7:
							G_num = (gx + 1) * grd->numGP2 + (gy + 1) * grd->numGP1 + (gz + 1);
							break;
						}
						wtr->O[i][dd].nearestgrdp = G_num;
						// store Oxygen coordinates
						wtr->O[i][dd].x[0] = x;
						wtr->O[i][dd].x[1] = y;
						wtr->O[i][dd].x[2] = z;

						// store Hydrogen coordinates SPC/TIP3P water
						fgets(line, lString, fi);
						//sscanf(line, "%*s%*d%*s%*s%*d%f%f%f", &x, &y, &z);
						GetChars(line, 30, 37, coor);
						x = atof(coor);
						GetChars(line, 38, 45, coor);
						y = atof(coor);
						GetChars(line, 46, 54, coor);
						z = atof(coor);

						wtr->H1[i][dd].x[0] = x;
						wtr->H1[i][dd].x[1] = y;
						wtr->H1[i][dd].x[2] = z;

						fgets(line, lString, fi);
						//sscanf(line, "%*s%*d%*s%*s%*d%f%f%f", &x, &y, &z);
						GetChars(line, 30, 37, coor);
						x = atof(coor);
						GetChars(line, 38, 45, coor);
						y = atof(coor);
						GetChars(line, 46, 54, coor);
						z = atof(coor);

						wtr->H2[i][dd].x[0] = x;
						wtr->H2[i][dd].x[1] = y;
						wtr->H2[i][dd].x[2] = z;

						if (par.site4_water) {
							fgets(line, lString, fi);
							//sscanf(line, "%*s%*d%*s%*s%*d%f%f%f", &x, &y, &z);
							GetChars(line, 30, 37, coor);
							x = atof(coor);
							GetChars(line, 38, 45, coor);
							y = atof(coor);
							GetChars(line, 46, 54, coor);
							z = atof(coor);

							wtr->H3[i][dd].x[0] = x;
							wtr->H3[i][dd].x[1] = y;
							wtr->H3[i][dd].x[2] = z;
						}

						if (par.site5_water) {
							fgets(line, lString, fi);
							//sscanf(line, "%*s%*d%*s%*s%*d%f%f%f", &x, &y, &z);
							GetChars(line, 30, 37, coor);
							x = atof(coor);
							GetChars(line, 38, 45, coor);
							y = atof(coor);
							GetChars(line, 46, 54, coor);
							z = atof(coor);

							wtr->H4[i][dd].x[0] = x;
							wtr->H4[i][dd].x[1] = y;
							wtr->H4[i][dd].x[2] = z;
						}

						//wtr->atn[i][dd] = atn;
						wtr->atn[i][dd] = fstWATatm + countAtm; //YY 4.21.2015
						//if (wtr->atn[i][dd] > 99999)//YY 4.21.2015
						// {
						// printf ("countAtm is %d\n", countAtm);
						// printf ("fstWATatm is %d\n", fstWATatm);
						// printf ("atom number now is %d = %d\n", wtr->atn[i][dd], fstWATatm + countAtm);
						// }
						wtr->resn[i][dd] = resn;
						dd++;
						//YY 12.6.2016 End of store water info inside binding site
					}
					if (par.site3_water) {
						countAtm += 3;
					}
					else if (par.site4_water) {
						countAtm += 4;
					}
					else if (par.site5_water) {
						countAtm += 5;
					}
				}

				if (strstr(line, "END"))
					break;
			}
			//printf("waters inside binding pocket: %d-%d\n", dd, cc);      
			wtr->WATinside[i] = dd;
		}

		wtr->largest_wat = 0;
		for (i = 0; i < numfr; i++) {
			if (wtr->WATinside[i] > wtr->largest_wat)
				wtr->largest_wat = wtr->WATinside[i];
		}

		fo = fopen("WATinside.mol2", "w");
		for (i = 0; i < numfr; i++) {
			fprintf(fo, "@<TRIPOS>MOLECULE\n****\n");
			fprintf(fo, "% 5d % 5d     0     0     0\n", wtr->WATinside[i] * 3,
					wtr->WATinside[i] * 2);
			fprintf(fo, "SMALL\nNO_CHARGES\n\n\n@<TRIPOS>ATOM\n");

			for (j = 0; j < wtr->WATinside[i]; j++) {
				fprintf(fo,
						"% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  1 <0>       % 8.4f\n",
						3 * j + 1, "WAT", wtr->O[i][j].x[0], wtr->O[i][j].x[1],
						wtr->O[i][j].x[2], "O", 0.0);
				fprintf(fo,
						"% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  1 <0>       % 8.4f\n",
						3 * j + 2, "WAT", wtr->H1[i][j].x[0],
						wtr->H1[i][j].x[1], wtr->H1[i][j].x[2], "H", 0.0);
				fprintf(fo,
						"% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  1 <0>       % 8.4f\n",
						3 * j + 3, "WAT", wtr->H2[i][j].x[0],
						wtr->H2[i][j].x[1], wtr->H2[i][j].x[2], "H", 0.0);
			}
			fprintf(fo, "@<TRIPOS>BOND\n");
			for (j = 0; j < wtr->WATinside[i]; j++) {
				fprintf(fo, "%6d %5d %5d %4s\n", 2 * j + 1, 3 * j + 1,
						3 * j + 2, "1");
				fprintf(fo, "%6d %5d %5d %4s\n", 2 * j + 2, 3 * j + 1,
						3 * j + 3, "1");
			}
			fprintf(fo, "\n\n");
		}
		fclose(fo);

		fo = fopen("WaterGrid.mol2", "w");
		fprintf(fo, "@<TRIPOS>MOLECULE\n****\n");
		fprintf(fo, "% 5d % 5d     0     0     0\n", grd->numGP3, 0);
		fprintf(fo, "SMALL\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n");

		denom = 1.0 / (sqrt_2 * sqrtPi * sigma * wtr->numfr);
		//printf("denom: %f\n", denom);
		for (num = 0; num < grd->numGP3; num++) {
			grd->WaterOccupancy[num] *= denom;      // normalize occupancy
			fprintf(fo,
					"% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  1 <0>       % 8.4f\n",
					num, "WAT", grd->GP_r[num][0], grd->GP_r[num][1],
					grd->GP_r[num][2], "O", grd->WaterOccupancy[num]);
		}
		fclose(fo);


		fo = fopen("grid_occupancy.dx", "w");
		fprintf(fo,
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
		for (num = 0; num < grd->numGP3; num++) {
			fprintf(fo, "%10.5f", grd->WaterOccupancy[num]);
			if ( (num+1)%3 == 0 ) {
				fprintf(fo, "\n");
			}
		}
		fclose(fo);

		// filter grid density > par.griddensity
		grid_counter = 0;
		for (gx = 0; gx < grd->NumGP_r[0]; gx++) {
			for (gy = 0; gy < grd->NumGP_r[1]; gy++) {
				for (gz = 0; gz < grd->NumGP_r[2]; gz++) {
					gnum = gx * grd->numGP2 + gy * grd->numGP1 + gz;
					if (grd->WaterOccupancy[gnum] > par.griddensity) {
						grid_counter++;
					}
				}
			}
		}
		fgrd->gridnum = grid_counter;

		fgrd->GP_r = (float **) calloc(fgrd->gridnum, sizeof(float *));
		fgrd->WaterOccupancy = (float*) calloc(fgrd->gridnum, sizeof(float));
		fgrd->ndx_in_grid = (int*) calloc(fgrd->gridnum, sizeof(float));
		for (i = 0; i < fgrd->gridnum; i++) {
			fgrd->GP_r[i] = (float *) calloc(3, sizeof(float));
			fgrd->WaterOccupancy[i] = 0;
			fgrd->ndx_in_grid[i] = 0;
		}

		new_counter = 0;
		for (gx = 0; gx < grd->NumGP_r[0]; gx++) {
			for (gy = 0; gy < grd->NumGP_r[1]; gy++) {
				for (gz = 0; gz < grd->NumGP_r[2]; gz++) {
					gnum = gx * grd->numGP2 + gy * grd->numGP1 + gz;
					if (grd->WaterOccupancy[gnum] > par.griddensity) {
						fgrd->GP_r[new_counter][0] = grd->GP_r[gnum][0];
						fgrd->GP_r[new_counter][1] = grd->GP_r[gnum][1];
						fgrd->GP_r[new_counter][2] = grd->GP_r[gnum][2];
						fgrd->WaterOccupancy[new_counter] =
								grd->WaterOccupancy[gnum];
						fgrd->ndx_in_grid[new_counter] = gnum;
						new_counter++;
					}
				}
			}
		}
		fo = fopen("FilterGrid.mol2", "w");
		fprintf(fo, "@<TRIPOS>MOLECULE\n****\n");
		fprintf(fo, "% 5d % 5d     0   0     0\n", fgrd->gridnum, 0);
		fprintf(fo, "SMALL\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n");

		for (num = 0; num < fgrd->gridnum; num++) {
			fprintf(fo,
					"% 7d %-4s    % 10.4f% 10.4f% 10.4f %-6s %4d WAT     % 8.4f\n",
					num, "WAT", fgrd->GP_r[num][0], fgrd->GP_r[num][1],
					fgrd->GP_r[num][2], "O", num, fgrd->WaterOccupancy[num]);
		}
		fclose(fo);
	}
}

//=================================================================================================
// Cluster grids using QT clustering
//=================================================================================================
void Cluster_QT(ligand *lig, grid *grd, water *wtr, hcluster *cluster_head) {
	QTClustering(grd, cluster_head, wtr);

	hcluster *cluster_0;
	cluster_0 = cluster_head->next;

	//OutputWaterIndex(cluster_0, wtr);
}

//=================================================================================================
// QT clustering algorithm
//=================================================================================================
void QTClustering(grid *grd, hcluster *cluster_head, water *wtr) {
	hcluster *cluster_0, *cluster, *prev_cluster;
	float size_water = par.size_water;
	float size_water2 = pow(size_water, 2);
	float sigma = size_water / 3.0;
	float sigma2 = pow(sigma, 2);
	float denom;
	float delta = grd->griddelta;       // grid spacing
	int del = (int) (size_water / delta) + 1; // number of gridpoints around center of water molecule up to
	int grid_cluster_dim = pow(2 * del, 3); //maximum nnumber of points could be included in each grid-centered cluster
	int cluster_found = 1;
	int gx, gy, gz, gnum;
	int k2, l2, m2, num;
	float x, y, z, xx, yy, zz, dd2;
	int mi, i;
	int largest_cluster_gnum;//grip point num of the cluster with the largest occupancy
	int num_largest_cluster_members;//number of grid points belong to the largest cluster
	float largest_cluster_occupancy;        //occupancy of the largest cluster
	float largest_cluster_avg_occupancy;
	int *largest_cluster_members;
	int *flag_in_cluster;       //flag of whether the grid has been clustered
	float *grid_cluster_occupancy;
	int **grid_cluster_members;
	int num_clusters;
	int mm, j, k, p, f, n, numwaters;
	float max_occu;
	float occupancy_cutoff = 8.0;
	int occupancy_cutoff2 = par.occupancy_cutoff;

	cluster_0 = NULL;
	flag_in_cluster = (int*) calloc(grd->numGP3, sizeof(int));
	grid_cluster_occupancy = (float*) calloc(grd->numGP3, sizeof(float));
	grid_cluster_members = (int**) calloc(grd->numGP3, sizeof(int*));
	for (i = 0; i < grd->numGP3; i++) {
		flag_in_cluster[i] = 0;
		grid_cluster_members[i] = (int*) calloc(grid_cluster_dim, sizeof(int));
	}
	largest_cluster_members = (int*) calloc(grid_cluster_dim, sizeof(int));

	num_clusters = 0;
	while (cluster_found) {
		largest_cluster_occupancy = -999.0;
		largest_cluster_avg_occupancy = -999.0;
		for (gx = 0; gx < grd->NumGP_r[0]; gx++) {
			for (gy = 0; gy < grd->NumGP_r[1]; gy++) {
				for (gz = 0; gz < grd->NumGP_r[2]; gz++) {
					gnum = gx * grd->numGP2 + gy * grd->numGP1 + gz;
					if (flag_in_cluster[gnum] == 0) {
						x = grd->GP_r[gnum][0];
						y = grd->GP_r[gnum][1];
						z = grd->GP_r[gnum][2];
						grid_cluster_occupancy[gnum] = 0.0;
						mi = 0;
						for (k2 = gx - del; k2 <= gx + del; k2++) {
							for (l2 = gy - del; l2 <= gy + del; l2++) {
								for (m2 = gz - del; m2 <= gz + del; m2++) {
									if (k2 >= 0 && k2 < grd->NumGP_r[0]
											&& l2 >= 0 && l2 < grd->NumGP_r[1]
											&& m2 >= 0
											&& m2 < grd->NumGP_r[2]) {
										num = k2 * grd->numGP2
												+ l2 * grd->numGP1 + m2;
										if (flag_in_cluster[num] == 0) {
											xx = grd->GP_r[num][0] - x;
											yy = grd->GP_r[num][1] - y;
											zz = grd->GP_r[num][2] - z;
											dd2 = xx * xx + yy * yy + zz * zz;
											if (dd2 <= size_water2) {
												grid_cluster_occupancy[gnum] +=
														grd->WaterOccupancy[num];
												grid_cluster_members[gnum][mi] =
														num;
												mi++;
											}
										}
									}
								}
							}
						}
						if (grid_cluster_occupancy[gnum]
								> largest_cluster_occupancy)
								//if(grid_cluster_occupancy[gnum]/mi > largest_cluster_avg_occupancy)
								{
							largest_cluster_gnum = gnum;
							largest_cluster_occupancy =
									grid_cluster_occupancy[gnum];
							largest_cluster_avg_occupancy =
									grid_cluster_occupancy[gnum] / mi;
							num_largest_cluster_members = mi;
							for (i = 0; i < num_largest_cluster_members; i++)
								largest_cluster_members[i] =
										grid_cluster_members[gnum][i];
						}
					}
				}
			}
		}
		printf("largest_cluster_occupancy: %f %f %d\n",
				largest_cluster_occupancy, largest_cluster_avg_occupancy,
				num_largest_cluster_members);
		//add into cluster list
		if (largest_cluster_occupancy > par.QT_occupancy_cutoff) {
			cluster = (hcluster *) calloc(1, sizeof(hcluster));
			if (cluster_0 == NULL) {
				cluster_0 = cluster;
			} else {
				prev_cluster->next = cluster;
			}
			cluster->numpoints = num_largest_cluster_members;
			cluster->center_occupancy =
					grd->WaterOccupancy[largest_cluster_gnum];
			cluster->overall_occupancy = largest_cluster_occupancy;
			cluster->center[0] = grd->GP_r[largest_cluster_gnum][0];
			cluster->center[1] = grd->GP_r[largest_cluster_gnum][1];
			cluster->center[2] = grd->GP_r[largest_cluster_gnum][2];
			cluster->grid_num = largest_cluster_gnum; //YY 1-15-2018 Forth number is the grid number in numGP3

			cluster->hydropoint = (float **) calloc(num_largest_cluster_members,
					sizeof(float *));
			cluster->hydroscore = (float *) calloc(num_largest_cluster_members,
					sizeof(float));
			cluster->grdpn = (int *) calloc(num_largest_cluster_members,
					sizeof(int));
			cluster->wtrindex = (int *) calloc(wtr->numfr, sizeof(int));

			for (j = 0; j < num_largest_cluster_members; j++) {
				cluster->hydropoint[j] = (float *) calloc(3, sizeof(float));
			}
			for (f = 0; f < wtr->numfr; f++) {
				cluster->wtrindex[f] = -1;
			}
			max_occu = -999.0;
			for (j = 0; j < num_largest_cluster_members; j++) {
				num = largest_cluster_members[j];
				for (p = 0; p < 3; p++) {
					cluster->hydropoint[j][p] = grd->GP_r[num][p];// atom numbers
				}
				cluster->hydroscore[j] = grd->WaterOccupancy[num];
				if (grd->WaterOccupancy[num] > max_occu)
					max_occu = grd->WaterOccupancy[num];
				cluster->grdpn[j] = num;

				//go over all water molecules in each frame, assign to the cluster
				for (f = 0; f < wtr->numfr; f++) {
					//nearestwtr = -1;
					for (n = 0; n < wtr->WATinside[f]; n++) {
						if (cluster->grdpn[j] == wtr->O[f][n].nearestgrdp) {
							cluster->wtrindex[f] = n;   //water index
							break;
						}
					}
				}
			}
			cluster->max_occupancy = max_occu;

			mm = 0;
			for (j = 0; j < wtr->numfr; j++) {
				if (cluster->wtrindex[j] >= 0) {
					mm++;
				}
			}
			cluster->wtrnum = mm;
			numwaters = mm;
			cluster->average_occupancy = largest_cluster_occupancy
					/ (num_largest_cluster_members);
			printf("max %f %f %f %d %d %d\n", max_occu,
					cluster->center_occupancy, cluster->average_occupancy,
					cluster->numpoints, num_largest_cluster_members, numwaters);

			prev_cluster = cluster;
			prev_cluster->next = NULL;

			//mark the grid point as clustered
			flag_in_cluster[largest_cluster_gnum] = 1;
			for (i = 0; i < num_largest_cluster_members; i++) {
				num = largest_cluster_members[i];
				flag_in_cluster[num] = 1;
			}
			num_clusters++;
		} else {
			cluster_found = 0;
		}
	}
	par.occupancy_cutoff = 200;

	cluster_0->numclusters = num_clusters;
	cluster_head->next = cluster_0; //give address of cluster_0 to the input pointer !!!this is very important!!!

	// free memory
	for (i = 0; i < grd->numGP3; i++) {
		free(grid_cluster_members[i]);
	}
	free(grid_cluster_members);
	free(flag_in_cluster);
	free(grid_cluster_occupancy);
	free(largest_cluster_members);

	//-----------------------OUTPUT cluster MOL2------------------------//
	FILE *fo;
	fo = fopen("HydrationSites.mol2", "w");
	fprintf(fo, "@<TRIPOS>MOLECULE\n****\n");
	fprintf(fo, "% 5d % 5d     0     0     0\n", num_clusters, 0);
	fprintf(fo, "SMALL min %f max %f\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n", 0.0,
			0.0);
	int ii = 0;
	for (cluster = cluster_0, i = 0; i < cluster_0->numclusters; cluster =
			cluster->next, i++) {
		fprintf(fo,
				"% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  %1d <0>       % 8.4f %5.2f\n",
				ii, "CLT", cluster->center[0], cluster->center[1],
				cluster->center[2], "O", i + 1, cluster->average_occupancy,
				(float) cluster->wtrnum / par.numfr);
		ii++;

	}
	fclose(fo);
}

//=================================================================================================
// OutputWaterIndex: output waterindex.wtndx file for calc enthalpy
//=================================================================================================
void OutputWaterIndex(hcluster *cluster_0, water *wtr) {
	int wtrln, frm, i, atn, add, resn, wtn, ci, cl, j;
	FILE *fo;
	int *atomlist, *resnlist, *cltnum, **cltlist, **frmlist;
	hcluster *cluster;
	int *cflag;
	printf("wtr->numfr %d*%d = %d\n", wtr->numfr, wtr->WATinside[0],
			wtr->WATinside[0] * wtr->numfr);

	atomlist = (int*) calloc(wtr->WATinside[0] * wtr->numfr, sizeof(int));
	resnlist = (int*) calloc(wtr->WATinside[0] * wtr->numfr, sizeof(int));
	cltnum = (int*) calloc(wtr->WATinside[0] * wtr->numfr, sizeof(int));

	cltlist = (int **) calloc(wtr->WATinside[0] * wtr->numfr, sizeof(int *));
	frmlist = (int **) calloc(wtr->WATinside[0] * wtr->numfr, sizeof(int *));

	cflag = (int *) calloc(cluster_0->numclusters, sizeof(int));
	int clusternum = cluster_0->numclusters;
	for (ci = 0; ci < clusternum; ci++) {
		cflag[ci] = 1;
	}

	wtrln = 0;
	for (cluster = cluster_0, cl = 0; cl < cluster_0->numclusters; cluster =
			cluster->next, cl++) {

		if (cflag[cl]) {
			//printf("cluster %d\n", cl);
			for (frm = 0; frm < wtr->numfr; frm++) {
				if (cluster->wtrindex[frm] >= 0) {
					cltlist[wtrln] = (int *) calloc(wtr->numfr, sizeof(int));
					frmlist[wtrln] = (int *) calloc(wtr->numfr, sizeof(int));

					wtn = cluster->wtrindex[frm];
					atn = wtr->atn[frm][wtn];
					resn = wtr->resn[frm][wtn];
					add = 1;
					for (i = 0; i < wtrln; i++) {
						if (atn == atomlist[i]) {
							add = 0;
							//if(cl == 33 || cl == 36) printf("clutnum %d %d %d %d\n", i, cltnum[i], cl, frm);
							cltlist[i][cltnum[i]] = cl;
							frmlist[i][cltnum[i]] = frm;
							cltnum[i]++;
							break;
						}
					}
					if (add) {
						atomlist[wtrln] = atn;
						resnlist[wtrln] = resn;
						//if(cl == 33 || cl == 36) printf("wtrln %d %d %d\n", wtrln, cl, frm);
						cltlist[wtrln][0] = cl;
						frmlist[wtrln][0] = frm;
						cltnum[wtrln] = 1;
						wtrln++;
					}
				}
			}
		}
	}
	fo = fopen("waterindex.wtndx", "w");
	fprintf(fo,
			"@wtr#  atom#   resn#   #of_clusters\n@             cluster#    frame#\n%6d %6d %6d\n",
			wtrln, cluster_0->numclusters, wtr->numfr);

	for (i = 0; i < wtrln; i++) {
		fprintf(fo, "%5d    %6d %5d   %4d\n", i, atomlist[i], resnlist[i],
				cltnum[i]);
		for (j = 0; j < cltnum[i]; j++) {
			fprintf(fo, "            %5d    %6d\n", cltlist[i][j], frmlist[i][j]);
		}
	}

	// free memory
	for (i = 0; i < wtr->WATinside[0] * wtr->numfr; i++) {
		free(frmlist[i]);
		free(cltlist[i]);
	}

	free(cltlist);
	free(frmlist);
	free(atomlist);

	free(resnlist);
	free(cltnum);
	fclose(fo);

}


//=================================================================================================
// Cluster grids using DBSCAN clustering
//=================================================================================================
void Cluster_DB(fgrid *fgrd, ligand *lig, grid *grd, water *wtr, hcluster *cluster_head) {
	dbScan(fgrd, grd, wtr, cluster_head);

	hcluster *cluster_0;
	cluster_0 = cluster_head->next;

	OutputWaterIndex(cluster_0, wtr);
}

//=================================================================================================
// DBSCAN clustering algorithm: output Cluster_DB.mol2 & HydrationSites.mol2
//=================================================================================================

void dbScan(fgrid *fgrd, grid *grd, water *wtr, hcluster *cluster_head) {

	printf("wtr->numfr %d*%d = %d\n", wtr->numfr, wtr->WATinside[0],
			wtr->WATinside[0] * wtr->numfr);

	int DATASET_SIZE = fgrd->gridnum;

	int FEATURES = 3; // dimensions of data

	float **data;
	int *clusters;
	int *clustered;
	int *visited;
	int *neigh_points;
	int i, j;

	data = (float**) calloc(DATASET_SIZE, sizeof(float*));
	for (i = 0; i < DATASET_SIZE; i++) {
		data[i] = (float*) calloc(FEATURES, sizeof(float));
	}

	for (i = 0; i < DATASET_SIZE; i++) {
		for (j = 0; j < FEATURES; j++) {
			data[i][j] = fgrd->GP_r[i][j];
		}
	}

	int next_cluster = 1, num_npoints;
	int MIN_POINTS = par.dbscan_START;
	clusters = (int*) calloc(DATASET_SIZE, sizeof(int));
	clustered = (int*) calloc(DATASET_SIZE, sizeof(int));

	while (MIN_POINTS >= par.dbscan_END) //CHANGE ME
	{
		printf("%d\n", MIN_POINTS);
		visited = (int*) calloc(DATASET_SIZE, sizeof(int));
		neigh_points = (int*) calloc(DATASET_SIZE * DATASET_SIZE, sizeof(int));

		for (i = 0; i < DATASET_SIZE; i++) {
			if (!visited[i] & !clustered[i]) {
				visited[i] = 1;
				num_npoints = regionQuery(0, i, DATASET_SIZE, FEATURES, data,
						neigh_points, MIN_POINTS);

				if (num_npoints > MIN_POINTS) {
					expandCluster(next_cluster, num_npoints, i, DATASET_SIZE,
							FEATURES, data, neigh_points, clusters, visited,
							MIN_POINTS, clustered);
					next_cluster++;
				}
			}
		}
		printf("%d\n", next_cluster);
		MIN_POINTS -= 10;
	}

	int num_clusters;
	int mm, qq, k, p, f, n, num, numwaters;
	float max_occu;
	float size_water = par.size_water;
	float size_water2 = pow(size_water, 2);
	float delta = grd->griddelta;       // grid spacing
	int del = (int) (size_water / delta) + 1; // number of gridpoints around center of water molecule up to
	int grid_cluster_dim = pow(2 * del, 3); //maximum nnumber of points could be included in each grid-centered cluster
	int *largest_cluster_members;
	int *WAT;

	int counter;
	float x, y, z, occu;
	float **save_cluster = (float**) calloc(next_cluster, sizeof(float*));
	for (i = 0; i < next_cluster; i++) {
		save_cluster[i] = (float*) calloc(5, sizeof(float));
	}

	hcluster *cluster_0, *prev_cluster, *cluster_new;
	cluster_0 = NULL;
	num_clusters = 0;

	for (j = 1; j < next_cluster; j++) {
		//printf("cluster index  j = %d\n", j);

		counter = 0;
		x = 0;
		y = 0;
		z = 0;
		occu = 0;
		largest_cluster_members = (int*) calloc(grid_cluster_dim, sizeof(int));

		for (i = 0; i < DATASET_SIZE; i++) {
			if (clusters[i] == j) {
				x += fgrd->GP_r[i][0];
				y += fgrd->GP_r[i][1];
				z += fgrd->GP_r[i][2];
				occu += fgrd->WaterOccupancy[i];
				//printf("ndx in grid   %d\n", fgrd->ndx_in_grid[i]);
				largest_cluster_members[counter] = fgrd->ndx_in_grid[i];
				counter++;
			}
		}
		//printf("number of gridpoints(%d) for counter j=%d\n", counter,j );

		//if( occu > par.QT_occupancy_cutoff)
		//{
		cluster_new = (hcluster *) calloc(1, sizeof(hcluster));
		if (cluster_0 == NULL) {
			cluster_0 = cluster_new;
		} else {
			prev_cluster->next = cluster_new;
		}

		cluster_new->numpoints = counter;

		save_cluster[num_clusters][4] = counter;
		save_cluster[num_clusters][0] = x / counter;
		save_cluster[num_clusters][1] = y / counter;
		save_cluster[num_clusters][2] = z / counter;
		save_cluster[num_clusters][3] = occu;

		cluster_new->overall_occupancy = occu;
		cluster_new->average_occupancy = occu / counter;
		cluster_new->center[0] = x / counter;
		cluster_new->center[1] = y / counter;
		cluster_new->center[2] = z / counter;

		//printf("occu %f, ave_occu %f, x %f, y %f, z %f\n", cluster_new->overall_occupancy, cluster_new->average_occupancy, cluster_new->center[0], cluster_new->center[1], cluster_new->center[2]);

		//printf("wtr->numfr = %d\n", wtr->numfr);
		WAT = (int*) calloc(counter, sizeof(int));

		cluster_new->wtrindex = (int *) calloc(wtr->numfr, sizeof(int));
		cluster_new->hydropoint = (float **) calloc(counter, sizeof(float *));
		cluster_new->hydroscore = (float *) calloc(counter, sizeof(float));
		cluster_new->grdpn = (int *) calloc(counter, sizeof(int));

		for (i = 0; i < counter; i++) {
			cluster_new->hydropoint[i] = (float *) calloc(3, sizeof(float));
		}
		//printf("iii = %d\n", i);
		for (f = 0; f < wtr->numfr; f++) {
			cluster_new->wtrindex[f] = -1;
		}

		max_occu = -999.0;
		//printf("max_occu %f\n", max_occu);

		printf("NOW cluster = %d\n", num_clusters);

		for (i = 0; i < counter; i++) {
			//printf("gridpoint index = %d\n", i);
			num = largest_cluster_members[i];
			//printf("gridpoint number = %d\n", num);

			for (p = 0; p < 3; p++) {
				//printf("xyz  %f\n", grd->GP_r[num][p]);
				cluster_new->hydropoint[i][p] = grd->GP_r[num][p];// atom numbers
			}
			cluster_new->hydroscore[i] = grd->WaterOccupancy[num];
			if (grd->WaterOccupancy[num] > max_occu) {
				max_occu = grd->WaterOccupancy[num];
				printf("%f\n", max_occu);
			}
			//printf("REACHED HERE\n");
			cluster_new->grdpn[i] = num;

			//go over all water molecules in each frame, assign to the cluster
			//mm = 0;
			for (f = 0; f < wtr->numfr; f++) {
				//nearestwtr = -1;
				mm = 0;

				for (n = 0; n < wtr->WATinside[f]; n++) {
					if (cluster_new->grdpn[i] == wtr->O[f][n].nearestgrdp) {
						cluster_new->wtrindex[f] = n;   //water index
						WAT[mm] = n;
						mm++;
						break;
					}
				}
				/*qq = mm;
				 printf("qq  %d\n", qq);
				 if (qq > 0)
				 {
				 printf("Cluster ndx %d   ", j-1);
				 printf("Grid ndx %d   ", i);
				 printf("frame %d   gnum = %d\n", f, num);

				 //printf("number of water near this gridpoint %d\n", mm);
				 //for (p=0; p<qq; p++)
				 //{
				 // printf("WAT[%d] index %d\n",p, WAT[p]);
				 //}

				 }
				 */
			}
		}
		cluster_new->max_occupancy = max_occu;

		mm = 0;
		for (f = 0; f < wtr->numfr; f++) {
			if (cluster_new->wtrindex[f] >= 0) {
				mm++;
			}
		}
		cluster_new->wtrnum = mm;
		numwaters = mm;
		printf("max %f %f %d %d %d\n", cluster_new->overall_occupancy,
				cluster_new->average_occupancy, cluster_new->numpoints, counter,
				numwaters);

		prev_cluster = cluster_new;
		prev_cluster->next = NULL;

		num_clusters++;
		//  }

	}

	printf("num_clusters %d\n", num_clusters);

	cluster_0->numclusters = num_clusters;
	cluster_head->next = cluster_0; //give address of cluster_0 to the input pointer !!!this is very important!!!
	// free memory
	free(largest_cluster_members);
	free(WAT);

	//-----------------------OUTPUT cluster MOL2------------------------//
	FILE *fo;
	fo = fopen("HydrationSites.mol2", "w");
	fprintf(fo, "@<TRIPOS>MOLECULE\n****\n");
	fprintf(fo, "% 5d % 5d     0     0     0\n", num_clusters, 0);
	fprintf(fo, "SMALL min %f max %f\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n", 0.0,
			0.0);
	int ii = 0;
	for (cluster_new = cluster_0, i = 0; i < cluster_0->numclusters;
			cluster_new = cluster_new->next, i++) {
		printf(
				"% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  %1d <0>       % 8.4f %5.2f\n",
				ii, "CLT", cluster_new->center[0], cluster_new->center[1],
				cluster_new->center[2], "O", i + 1,
				cluster_new->overall_occupancy,
				(float) cluster_new->wtrnum / par.numfr);
		fprintf(fo,
				"% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  %1d <0>       % 8.4f %5.2f\n",
				ii, "CLT", cluster_new->center[0], cluster_new->center[1],
				cluster_new->center[2], "O", i + 1,
				cluster_new->overall_occupancy,
				(float) cluster_new->wtrnum / par.numfr);
		ii++;

	}
	fclose(fo);

	//-----------------------OUTPUT cluster MOL2------------------------//

	/*FILE  *fo;
	 fo = fopen( "Cluster_DBSCAN.mol2", "w");
	 fprintf(fo, "@<TRIPOS>MOLECULE\n****\n");
	 fprintf(fo, "% 5d % 5d     0   0     0\n", next_cluster-1, 0);
	 fprintf(fo, "SMALL\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n");
	 for (j=0; j<num_clusters; j++)
	 {
	 fprintf(fo, "% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  %1d <0>     % 8.4f\n", j, "WAT", save_cluster[j][0], save_cluster[j][1], save_cluster[j][2], "O", j+1,save_cluster[j][3]);
	 }
	 fclose(fo);


	 mkdir("DBSCAN_mol2", S_IRWXU | S_IRWXG | S_IRWXO);
	 chdir("DBSCAN_mol2");
	 for (j=0; j<next_cluster-1; j++)
	 {
	 char fname[32];
	 snprintf(fname, sizeof(char) * 32, "Cluster_DBSCAN_%i.mol2", j);
	 fo = fopen( fname, "w");
	 fprintf(fo, "@<TRIPOS>MOLECULE\n****\n");
	 fprintf(fo, "% 5d % 5d     0   0     0\n", (int)save_cluster[j][4], 0);
	 fprintf(fo, "SMALL\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n");

	 int counter = 0;
	 for(i = 0; i < DATASET_SIZE; i++)
	 {
	 if (clusters[i] == j+1)
	 {
	 fprintf(fo, "% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  1 <0>     % 8.4f\n", counter, "WAT", fgrd->GP_r[i][0], fgrd->GP_r[i][1], fgrd->GP_r[i][2], "O", fgrd->WaterOccupancy[i]);
	 counter ++;
	 }
	 }
	 // printf("%f  %f  %f\n", save_cluster[j][0]/counter, save_cluster[j][1]/counter, save_cluster[j][2]/counter);
	 fclose(fo);

	 }
	 chdir("../");
	 */

}

void expandCluster(float cluster_no, int num_npoints, int index,
		int DATASET_SIZE, int FEATURES, float **data, int *neigh_points,
		int *clusters, int *visited, int MIN_POINTS, int *clustered) {
	clusters[index] = cluster_no;
	clustered[index]++;
	data[index][0] = 0;
	data[index][1] = 0;
	data[index][2] = 0;

	int i, count = 0;

	for (i = 0; i < num_npoints; i++) {
		if (!visited[neigh_points[i]]) {
			visited[neigh_points[i]] = 1;

			count = regionQuery(num_npoints, neigh_points[i], DATASET_SIZE,
					FEATURES, data, neigh_points, MIN_POINTS);

			if (count >= MIN_POINTS) {
				num_npoints += count;
			}
		}

		if (!clusters[neigh_points[i]]) {
			clusters[neigh_points[i]] = cluster_no;
			data[neigh_points[i]][0] = 0;
			data[neigh_points[i]][1] = 0;
			data[neigh_points[i]][2] = 0;
			clustered[neigh_points[i]]++;
		}
	}
}

int regionQuery(int start, int index, int DATASET_SIZE, int FEATURES,
		float **data, int *neigh_points, int MIN_POINTS) {
	int i, j, count = 0;
	float distance, temp;

	//printf("function regionQuery\n");
	for (i = 0; i < DATASET_SIZE; i++) {
		if (i != index) {
			distance = 0;

			for (j = 0; j < FEATURES; j++) {
				temp = data[i][j] - data[index][j];

				distance += temp * temp;
			}

			if (distance <= EPSILON) {
				//printf("distance %f\n", distance );
				neigh_points[start + count] = i;
				count++;
			}
		}
	}
	return count;
}



void Grid_cluster(grid *grd, water *wtr, hcluster *cluster_head) {
	hcluster *cluster_0, *cluster, *prev_cluster;
	float size_water = par.size_water;
	float size_water2 = pow(size_water, 2);
	float sigma = size_water / 3.0;
	float sigma2 = pow(sigma, 2);
	float denom;
	float delta = grd->griddelta;       // grid spacing
	int del = (int) (size_water / delta) + 1; // number of grid points around center of water molecule up to
	int grid_cluster_dim = pow(2 * del, 3); //maximum number of points could be included in each grid-centered cluster


	float cluster_bulk[3] = { 0.1, 0.5, 0.9 };
	int numGridSplit = 3;

	int gx, gy, gz, gnum;
	int k2, l2, m2, num;
	float x, y, z, xx, yy, zz, dd2;
	int mi, i, ii, jj, kk;

	int *largest_cluster_members;
	float *grid_cluster_occupancy;
	int **grid_cluster_members;

	int mm, j, k, p, f, n;
	float max_occu;

	cluster_0 = NULL;
	//grid_cluster_occupancy = (float*) calloc(grd->numGP3, sizeof(float));
	grid_cluster_members = (int**) calloc(grd->numGP3, sizeof(int*));
	for (i = 0; i < grd->numGP3; i++) {
		grid_cluster_members[i] = (int*) calloc(grid_cluster_dim, sizeof(int));
	}
	//largest_cluster_members = (int*) calloc(grid_cluster_dim, sizeof(int));

	int num_clusters = 0;
	int zeroGrid = 0;

	for (gx = 0; gx < grd->NumGP_r[0]; gx++) {
		for (gy = 0; gy < grd->NumGP_r[1]; gy++) {
			for (gz = 0; gz < grd->NumGP_r[2]; gz++) {
				gnum = gx * grd->numGP2 + gy * grd->numGP1 + gz;
				if ( grd->WaterOccupancy[gnum] > 0.045 ) {
					x = grd->GP_r[gnum][0];
					y = grd->GP_r[gnum][1];
					z = grd->GP_r[gnum][2];

					//grid_cluster_occupancy[gnum] += grd->WaterOccupancy[gnum];

					mi = 0;
					for (k2 = gx - del; k2 <= gx + del; k2++) {
						for (l2 = gy - del; l2 <= gy + del; l2++) {
							for (m2 = gz - del; m2 <= gz + del; m2++) {
								if (k2 >= 0 && k2 < grd->NumGP_r[0]
								 && l2 >= 0 && l2 < grd->NumGP_r[1]
								 && m2 >= 0 && m2 < grd->NumGP_r[2])
								{
									num = k2 * grd->numGP2 + l2 * grd->numGP1 + m2;
									xx = grd->GP_r[num][0] - x;
									yy = grd->GP_r[num][1] - y;
									zz = grd->GP_r[num][2] - z;
									dd2 = xx * xx + yy * yy + zz * zz;
									if (dd2 <= size_water2) {
										//grid_cluster_occupancy[gnum] += grd->WaterOccupancy[num];
										grid_cluster_members[gnum][mi] = num;
										mi++;
									}
								}
							}
						}
					}
					//for (i = 0; i < mi; i++)
					//    largest_cluster_members[i] = grid_cluster_members[gnum][i];
					//printf("largest_cluster_members: %d\n", mi);

					//add into cluster list
					cluster = (hcluster *) calloc(1, sizeof(hcluster));
					if (cluster_0 == NULL) {
						cluster_0 = cluster;
					} else {
						prev_cluster->next = cluster;
					}
					cluster->numpoints = mi;
					cluster->center_occupancy = grd->WaterOccupancy[gnum];
					cluster->center[0] = grd->GP_r[gnum][0];
					cluster->center[1] = grd->GP_r[gnum][1];
					cluster->center[2] = grd->GP_r[gnum][2];
					cluster->grid_num = gnum; //YY 1-15-2018 Forth number is the grid number in numGP3
					//printf("grid num: %d %d\n", cluster->grid_num, gnum);

///*
					cluster->wtrindex = (int *) calloc(wtr->numfr, sizeof(int));
					for (f = 0; f < wtr->numfr; f++) {
						cluster->wtrindex[f] = -1;
					}

					// if grid center is the nearest grid of a water molecule... 
					for (f = 0; f < wtr->numfr; f++) {
						bool stop = false;
						for (n = 0; n < wtr->WATinside[f] && !stop; n++) {
							if (gnum == wtr->O[f][n].nearestgrdp) {
								cluster->wtrindex[f] = n;   //water index
								//printf("gnum=%d  frame=%d  wtrindex=%d \n", gnum, f, n);
								stop = true;
							}

							for (j = 0; (j < mi) && !stop; j++) {
								num = grid_cluster_members[gnum][j];
								if (num == wtr->O[f][n].nearestgrdp) {
									cluster->wtrindex[f] = n;   //water index
									//printf("   num=%d  frame=%d  wtrindex=%d \n", num, f, n);
									stop = true;
								}
							}
						}
					}

					mm = 0;
					for (j = 0; j < wtr->numfr; j++) {
						if (cluster->wtrindex[j] >= 0) {
							mm++;
						}
					}
					cluster->wtrnum = mm;
					prev_cluster = cluster;
					prev_cluster->next = NULL;

					num_clusters++;
//*/                    
/*
					// YY 3-23-18 rewrite only useful info to speed up grid calculation 

					cluster->hydropoint = (float **) calloc( mi, sizeof(float *));
					cluster->hydroscore = (float *) calloc( mi, sizeof(float));
					cluster->grdpn = (int *) calloc(mi, sizeof(int));
					cluster->wtrindex = (int *) calloc(wtr->numfr, sizeof(int));

					for (j = 0; j < mi; j++) {
						cluster->hydropoint[j] = (float *) calloc(3, sizeof(float));
					}
					for (f = 0; f < wtr->numfr; f++) {
						cluster->wtrindex[f] = -1;
					}
					max_occu = -999.0;
					for (j = 0; j < mi; j++) {
						num = grid_cluster_members[gnum][j];
						for (p = 0; p < 3; p++) {
							cluster->hydropoint[j][p] = grd->GP_r[num][p];// atom numbers
						}
						cluster->hydroscore[j] = grd->WaterOccupancy[num];
						if (grd->WaterOccupancy[num] > max_occu)
							max_occu = grd->WaterOccupancy[num];
						cluster->grdpn[j] = num;

						//go over all water molecules in each frame, assign to the cluster
						for (f = 0; f < wtr->numfr; f++) {
							//nearestwtr = -1;
							for (n = 0; n < wtr->WATinside[f]; n++) {
								if (cluster->grdpn[j] == wtr->O[f][n].nearestgrdp) {
									cluster->wtrindex[f] = n;   //water index
									break;
								}
							}
						}
					}
					cluster->max_occupancy = max_occu;
					mm = 0;
					for (j = 0; j < wtr->numfr; j++) {
						if (cluster->wtrindex[j] >= 0) {
							mm++;
						}
					}
					cluster->wtrnum = mm;
					//printf("max_occu %f, cen_occu %f, cluster num grids %d, cluster water number %d\n",
					//        max_occu, cluster->center_occupancy, cluster->numpoints, cluster->wtrnum );
					prev_cluster = cluster;
					prev_cluster->next = NULL;

					num_clusters++;
					//printf("grid calulation for %d\r", num_clusters);
*/
				}
				else{
					zeroGrid ++; 
				}
				printf("Grid calculated %d -- Zero grid %d \r", num_clusters, zeroGrid );
			}
		}
	}

	cluster_0->numclusters = num_clusters;
	cluster_head->next = cluster_0; //give address of cluster_0 to the input pointer !!!this is very important!!!

	// free memory
	for (i = 0; i < grd->numGP3; i++) {
		free(grid_cluster_members[i]);
	}
	free(grid_cluster_members);
	//free(grid_cluster_occupancy);
	//free(largest_cluster_members);

}


//=================================================================================================
// Bulk of DBSCAN Clutering
//=================================================================================================
void BulkCluster(fgrid *fgrd, grid *grd, water *wtr, hcluster *cluster_head) {
	hcluster *cluster_0, *cluster, *prev_cluster;
	float size_water = par.size_water;
	float size_water2 = pow(size_water, 2);
	float sigma = size_water / 3.0;
	float sigma2 = pow(sigma, 2);
	float denom;
	float delta = grd->griddelta;       // grid spacing
	int del = (int) (size_water / delta) + 1; // number of gridpoints around center of water molecule up to
	int grid_cluster_dim = pow(2 * del, 3); //maximum nnumber of points could be included in each grid-centered cluster

	float cluster_found;
	//float cluster_bulk[5] = {0.25, 0.333, 0.5, 0.667, 0.75};
	float cluster_bulk[3] = { 0.1, 0.5, 0.9 };
	int numGridSplit = 3;

	int gx, gy, gz, gnum;
	int k2, l2, m2, num;
	float x, y, z, xx, yy, zz, dd2;
	int mi, i, ii, jj, kk;
	int largest_cluster_gnum;//grip point num of the cluster with the largest occupancy
	int num_largest_cluster_members;//number of grid points belong to the largest cluster
	float largest_cluster_occupancy;        //occupancy of the largest cluster
	float largest_cluster_avg_occupancy;
	int *largest_cluster_members;
	int *flag_in_cluster;       //flag of whether the grid has been clustered
	float *grid_cluster_occupancy;
	int **grid_cluster_members;
	int num_clusters;
	int mm, j, k, p, f, n, numwaters;
	float max_occu;
	float occupancy_cutoff = 8.0;
	int occupancy_cutoff2 = par.occupancy_cutoff;

	cluster_0 = NULL;
	flag_in_cluster = (int*) calloc(grd->numGP3, sizeof(int));
	grid_cluster_occupancy = (float*) calloc(grd->numGP3, sizeof(float));
	grid_cluster_members = (int**) calloc(grd->numGP3, sizeof(int*));
	for (i = 0; i < grd->numGP3; i++) {
		flag_in_cluster[i] = 0;
		grid_cluster_members[i] = (int*) calloc(grid_cluster_dim, sizeof(int));
	}
	largest_cluster_members = (int*) calloc(grid_cluster_dim, sizeof(int));

	num_clusters = 0;

	for (ii = 0; ii < numGridSplit; ii++) {
		for (jj = 0; jj < numGridSplit; jj++) {
			for (kk = 0; kk < numGridSplit; kk++) {

				//cluster_found = cluster_bulk[ii];
				//printf("Cluster_found = %f\n", cluster_found);

				largest_cluster_occupancy = -999.0;
				largest_cluster_avg_occupancy = -999.0;

				//float yy1 = cluster_found + 0.1;
				//float yy2 = cluster_found - 0.1;
				gx = round(grd->NumGP_r[0] * cluster_bulk[ii]);
				gy = round(grd->NumGP_r[1] * cluster_bulk[jj]);
				gz = round(grd->NumGP_r[2] * cluster_bulk[kk]);
				gnum = gx * grd->numGP2 + gy * grd->numGP1 + gz;
				printf("gx gy gz gnum %d %d %d %d\n", gx, gy, gz, gnum);

				x = grd->GP_r[gnum][0];
				y = grd->GP_r[gnum][1];
				z = grd->GP_r[gnum][2];
				grid_cluster_occupancy[gnum] = 0.0;
				mi = 0;
				for (k2 = gx - del; k2 <= gx + del; k2++) {
					for (l2 = gy - del; l2 <= gy + del; l2++) {
						for (m2 = gz - del; m2 <= gz + del; m2++) {
							if (k2 >= 0 && k2 < grd->NumGP_r[0] && l2 >= 0
									&& l2 < grd->NumGP_r[1] && m2 >= 0
									&& m2 < grd->NumGP_r[2]) {
								num = k2 * grd->numGP2 + l2 * grd->numGP1 + m2;
								if (flag_in_cluster[num] == 0) {
									xx = grd->GP_r[num][0] - x;
									yy = grd->GP_r[num][1] - y;
									zz = grd->GP_r[num][2] - z;
									dd2 = xx * xx + yy * yy + zz * zz;
									if (dd2 <= size_water2) {
										grid_cluster_occupancy[gnum] +=
												grd->WaterOccupancy[num];
										grid_cluster_members[gnum][mi] = num;
										mi++;
									}
								}
							}
						}
					}
				}
				if (grid_cluster_occupancy[gnum] > largest_cluster_occupancy)
				//if(grid_cluster_occupancy[gnum]/mi > largest_cluster_avg_occupancy)
						{
					largest_cluster_gnum = gnum;
					largest_cluster_occupancy = grid_cluster_occupancy[gnum];
					largest_cluster_avg_occupancy = grid_cluster_occupancy[gnum]
							/ mi;
					num_largest_cluster_members = mi;
					for (i = 0; i < num_largest_cluster_members; i++)
						largest_cluster_members[i] =
								grid_cluster_members[gnum][i];
				}

				printf("largest_cluster_occupancy: %f %f %d\n",
						largest_cluster_occupancy,
						largest_cluster_avg_occupancy,
						num_largest_cluster_members);

				//add into cluster list
				if (largest_cluster_occupancy > 0) {
					cluster = (hcluster *) calloc(1, sizeof(hcluster));
					if (cluster_0 == NULL) {
						cluster_0 = cluster;
					} else {
						prev_cluster->next = cluster;
					}
					cluster->numpoints = num_largest_cluster_members;
					cluster->center_occupancy =
							grd->WaterOccupancy[largest_cluster_gnum];
					cluster->overall_occupancy = largest_cluster_occupancy;
					cluster->center[0] = grd->GP_r[largest_cluster_gnum][0];
					cluster->center[1] = grd->GP_r[largest_cluster_gnum][1];
					cluster->center[2] = grd->GP_r[largest_cluster_gnum][2];

					cluster->hydropoint = (float **) calloc(
							num_largest_cluster_members, sizeof(float *));
					cluster->hydroscore = (float *) calloc(
							num_largest_cluster_members, sizeof(float));
					cluster->grdpn = (int *) calloc(num_largest_cluster_members,
							sizeof(int));
					cluster->wtrindex = (int *) calloc(wtr->numfr, sizeof(int));

					for (j = 0; j < num_largest_cluster_members; j++) {
						cluster->hydropoint[j] = (float *) calloc(3,
								sizeof(float));
					}
					for (f = 0; f < wtr->numfr; f++) {
						cluster->wtrindex[f] = -1;
					}
					max_occu = -999.0;
					for (j = 0; j < num_largest_cluster_members; j++) {
						num = largest_cluster_members[j];
						for (p = 0; p < 3; p++) {
							cluster->hydropoint[j][p] = grd->GP_r[num][p];// atom numbers
						}
						cluster->hydroscore[j] = grd->WaterOccupancy[num];
						if (grd->WaterOccupancy[num] > max_occu)
							max_occu = grd->WaterOccupancy[num];
						cluster->grdpn[j] = num;

						//go over all water molecules in each frame, assign to the cluster
						for (f = 0; f < wtr->numfr; f++) {
							//nearestwtr = -1;
							for (n = 0; n < wtr->WATinside[f]; n++) {
								if (cluster->grdpn[j]
										== wtr->O[f][n].nearestgrdp) {
									cluster->wtrindex[f] = n;   //water index
									break;
								}
							}
						}
					}
					cluster->max_occupancy = max_occu;

					mm = 0;
					for (j = 0; j < wtr->numfr; j++) {
						if (cluster->wtrindex[j] >= 0) {
							mm++;
						}
					}
					cluster->wtrnum = mm;
					numwaters = mm;
					cluster->average_occupancy = largest_cluster_occupancy
							/ (num_largest_cluster_members);
					printf("max %f %f %f %d %d %d\n", max_occu,
							cluster->center_occupancy,
							cluster->average_occupancy, cluster->numpoints,
							num_largest_cluster_members, numwaters);

					prev_cluster = cluster;
					prev_cluster->next = NULL;

					//mark the grid point as clustered
					flag_in_cluster[largest_cluster_gnum] = 1;
					for (i = 0; i < num_largest_cluster_members; i++) {
						num = largest_cluster_members[i];
						flag_in_cluster[num] = 1;
					}
					num_clusters++;
				}
				printf("num_clusters %d\n", num_clusters);

			}
		}
	}
	par.occupancy_cutoff = 200;

	cluster_0->numclusters = num_clusters;
	cluster_head->next = cluster_0; //give address of cluster_0 to the input pointer !!!this is very important!!!

	// free memory
	for (i = 0; i < grd->numGP3; i++) {
		free(grid_cluster_members[i]);
	}
	free(grid_cluster_members);
	free(flag_in_cluster);
	free(grid_cluster_occupancy);
	free(largest_cluster_members);

	//-----------------------OUTPUT cluster MOL2------------------------//
	FILE *fo;
	fo = fopen("HydrationSites.mol2", "w");
	fprintf(fo, "@<TRIPOS>MOLECULE\n****\n");
	fprintf(fo, "% 5d % 5d     0     0     0\n", num_clusters, 0);
	fprintf(fo, "SMALL min %f max %f\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n", 0.0,
			0.0);
	ii = 0;
	for (cluster = cluster_0, i = 0; i < cluster_0->numclusters; cluster =
			cluster->next, i++) {
		fprintf(fo,
				"% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  %1d <0>       % 8.4f %5.2f\n",
				ii, "CLT", cluster->center[0], cluster->center[1],
				cluster->center[2], "O", i + 1, cluster->average_occupancy,
				(float) cluster->wtrnum / par.numfr);
		ii++;

	}
	fclose(fo);
}

