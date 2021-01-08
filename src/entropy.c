/***************************************************************************
 *   This script calculates water entropy based on water occupancy clusters*
 *   1) calculates Euler Angles                                            *
 *   2) computes 6x6 (translation+rotation) or 3x3 covariance matrix for   *
 *   PCA.                                                                  *
 *   3) computes probablity density function (70 bins)                     *
 *   4) does simpsons' integral on each dimension                          *
 ***************************************************************************/

#include "entropy.h"

//=================================================================================================
// EulerAngle: generate Euler angles for wates
//=================================================================================================
void CalcEntropy(water *wtr, hcluster *cluster_0, energy *engy)
{
	//calculate Euler angles for each water inside cluster
	printf("\ncalculating Euler Angle for rotation...\n");
	EulerAngle(wtr, cluster_0);

	//calculate 6*6 covariance matrix
	if(par.covdimension == 6 || par.covdimension == 9)
	{
		mkdir("6x6", S_IRWXU | S_IRWXG | S_IRWXO);
		chdir("6x6");

		CovarianceMatrix6x6(wtr, cluster_0, engy);
		
		chdir("../");
	}
	//calculate 3x3 covariance matrix 
	if(par.covdimension == 3 || par.covdimension == 9)
	{
		mkdir("3x3", S_IRWXU | S_IRWXG | S_IRWXO);
		chdir("3x3");

		CovarianceMatrix3x3(wtr, cluster_0, engy);

		chdir("../");
	}
}

//=================================================================================================
// EulerAngle: generate Euler angles for wates
//=================================================================================================
void EulerAngle(water *wtr, hcluster *cluster_0)
{
	//go over each cluster, calculate euler angle for waters inside of the cluster 
	//store euler angle under cluster 
	hcluster	*cluster;
	int			i, frm, wtn, j;
	float		v1[3], v2[3], v3[3], tmp[3], ss;
	float		alpha, beta, gamma, sa, ca, sb, cb, sg, cg, v1dotv2;
	
	for(cluster = cluster_0, j = 0; j < cluster_0->numclusters; cluster = cluster->next, j++)
	{
		cluster->wtrEulerAgl = (float **)calloc(wtr->numfr, sizeof(float*));
		for(frm = 0; frm < wtr->numfr; frm++)
		{
			cluster->wtrEulerAgl[frm] = (float *)calloc(3, sizeof(float));
		}
		if(par.grid_energy == 0) printf("cluster %d: average_occupancy: %f\n", j, cluster->average_occupancy);
		// go over each water in the cluster
		for(frm = 0; frm < wtr->numfr; frm++)
		{
			if(cluster->wtrindex[frm] >= 0)
			{
				wtn = cluster->wtrindex[frm];
				// create normal vectors
				
				CreateNormVector(wtr->H1[frm][wtn].x, wtr->O[frm][wtn].x, v1);
				
				CreateNormVector(wtr->H2[frm][wtn].x, wtr->O[frm][wtn].x, tmp);
				
				VdotV(v1, tmp, &v1dotv2);
				for(i = 0; i < 3; i++)
				{
					//v2[i] = tmp[i] + sinf(14.52*Pi/180.0)*v1[i];
					v2[i] = tmp[i] - v1dotv2*v1[i];	//gram-schmidt diagonilization
				}
				
				NormalizeVector(v2);
				
				VxV(v1, v2, v3);

/*				VdotV(v1, v2, &ss);
				printf("test v1 v2 ss %f\n", ss);
				VdotV(v1, v3, &ss);
				printf("test v1 v3 ss %f\n", ss);
				VdotV(v2, v3, &ss);
				printf("test v2 v3 ss %f\n", ss);
*/
				// calculate Euler Angles 
				sb = sqrtf(v3[0]*v3[0]+v3[1]*v3[1]);
				cb = v3[2];
				beta = atan2f(sb, cb);
				alpha = atan2f(v3[0], -v3[1]);
				gamma = atan2f(v1[2], v2[2]);
				
				// store Euler angles 
				cluster->wtrEulerAgl[frm][0] = alpha;
				cluster->wtrEulerAgl[frm][1] = beta;
				cluster->wtrEulerAgl[frm][2] = gamma;
			}
		}
		
	}
}

//====================================================================================================
// CovarianceMatrix: generate zero-mean covariance matrix for scaled positions and Euler angles (6x6)
//====================================================================================================
void CovarianceMatrix6x6(water *wtr, hcluster *cluster_0, energy *engy)
{
	hcluster	*cluster;
	float		**coor, **cov, coor_mean[6], **newcoor, newcov[6][6];
	int			m, n, s, i, j, frm, wtn, cl;
	float		egvl[6], e[6];
	FILE		*fo;
	
	printf("6x6 covariance matrix...\n");
	//go over each cluster 
	cov = (float **)calloc(6, sizeof(float*));
	for(i = 0; i < 6; i++)
	{
		cov[i] = (float*)calloc(6, sizeof(float));
	}
	for(cluster = cluster_0, cl = 0; cl < cluster_0->numclusters; cluster = cluster->next, cl++)
	{
		//assign coordinate and euler angle to coor
		coor = (float**)calloc(cluster->wtrnum, sizeof(float*));
		for(i = 0; i < cluster->wtrnum; i++)
		{
			coor[i] = (float*)calloc(6, sizeof(float));
		}
		for(i = 0; i < 6; i++)
		{
			coor_mean[i] = 0.0;
		}
		for(m = 0; m < 6; m++)
		{
			for(n = 0; n < 6; n++)
			{
				cov[m][n] = 0.0;
				newcov[m][n]=0.0;
			}
		}

		m = 0;
		for(frm = 0; frm < wtr->numfr; frm++)
		{
			if(cluster->wtrindex[frm] >= 0)
			{
				wtn = cluster->wtrindex[frm];
				for(i = 0; i < 3; i++)
				{
					coor[m][i] = wtr->O[frm][wtn].x[i];
					coor[m][i+3] = cluster->wtrEulerAgl[frm][i];
				}
				m++;
			}
		}
		if(par.grid_energy == 0) printf("water inside cluster%d: %d (=%d)\n", cl, m, cluster->wtrnum);
		
		//calculate mean 
		for(m = 0; m < cluster->wtrnum; m++)
		{	
			for(i = 0; i < 6; i++)
			{
				coor_mean[i] += coor[m][i];
			}
		}		
		
		for(i = 0; i < 6; i++)
		{
			coor_mean[i] /= cluster->wtrnum;
		}

		for(i = 0; i < cluster->wtrnum; i++)
		{
			for(m = 0; m < 6; m++)
			{
				coor[i][m] -= coor_mean[m];
			}				
		}

		//construct covariance matrix
		for(m = 0; m < 6; m++)
		{
			for(n = 0; n < 6; n++)
			{
				for(i = 0; i < cluster->wtrnum; i++)
				{
					cov[m][n] += coor[i][m]*coor[i][n];
				}				
			}
		}

		//printf("covariance matrix:\n");
		for(m = 0; m < 6; m++)
		{
			for(n = 0; n < 6; n++)
			{
				if(cluster->wtrnum == 1)	cov[m][n] /= cluster->wtrnum;
				else 						cov[m][n] /= (cluster->wtrnum - 1);  // ???? n or n-1???
				cluster->cov[m][n] = cov[m][n];
			}
			//printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f \n", cluster->cov[m][0], cluster->cov[m][1], cluster->cov[m][2], cluster->cov[m][3], cluster->cov[m][4], cluster->cov[m][5] );
		}	
		
		//calculate eigenvector and eigenvalue
		tred2(cov, 6, egvl, e);

		tqli(egvl, e, 6, cov);

		//printf("eigenvectors:\n");
		for(m = 0; m < 6; m++)
		{
			for(n = 0; n < 6; n++)
			{
				cluster->eigenvt[m][n] = cov[m][n];
			}
			cluster->eigenvl[m] = egvl[m];
			//printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f: %f\n", cov[0][m], cov[1][m], cov[2][m], cov[3][m], cov[4][m], cov[5][m], egvl[m]);
		}

		newcoor = (float**)calloc(6, sizeof(float*));
		for(m = 0; m < 6; m++)
		{
			newcoor[m] = (float*)calloc(cluster->wtrnum, sizeof(float));
			for(n = 0; n < cluster->wtrnum; n++)
				newcoor[m][n] = 0.0;
		}
		
		for(s = 0; s < cluster->wtrnum; s++)
		{
			for(m = 0; m < 6; m++)
			{
				for(n = 0; n < 6; n++)
				{
					if(n == 4)	//Euler angle beta
					{
						newcoor[m][s] += cov[n][m] * coor[s][n];// * sin(coor[s][n]);
					}
					else
					{
						newcoor[m][s] += cov[n][m] * coor[s][n];
					}
				}
			}
		}

		cluster->PCAcoor = (float**)calloc(6, sizeof(float*));
		for(m = 0; m < 6; m++)
		{
			cluster->PCAcoor[m] = (float*)calloc(cluster->wtrnum, sizeof(float));
			for(s = 0; s < cluster->wtrnum; s++)
			{
				cluster->PCAcoor[m][s] = newcoor[m][s];
			}
		}
		
		//verify calculation by computing new covariance 
		//if it is correct, new covariance matrix should be diagonal (with eigenvalues in the diagonal entries)
/*		for(m = 0; m < 6; m++)
		{
			for(n = 0; n < 6; n++)
			{
				for(i = 0; i < cluster->wtrnum; i++)
				{
					//cov[m][n] += (coor[i][m] - coor_mean[m])*(coor[i][n] - coor_mean[n]);
					newcov[m][n] += newcoor[m][i]*newcoor[n][i];
				}				
			}
		}
		printf("new covariance matrix:\n");
		for(m = 0; m < 6; m++)
		{
			for(n = 0; n < 6; n++)
			{
				newcov[m][n] /= (cluster->wtrnum - 1);  // ???? n or n-1???
			}
			printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f \n", newcov[m][0], newcov[m][1], newcov[m][2], newcov[m][3], newcov[m][4], newcov[m][5] );
		}	
*/
		//free memory
		for(i = 0; i < cluster->wtrnum; i++)
		{
			free(coor[i]);
		}
		free(coor);	
		for(i = 0; i < 6; i++)
		{
			free(newcoor[i]);
		}
		free(newcoor);
	}
	
	for(m = 0; m < 6; m++)
		free(cov[m]);
	free(cov);
	
	//-----------------OUTPUT file for histogram---------------//
	fo = fopen("PCAhistogram6x6.txt", "w");
	fprintf(fo, "This file has the coordinates projected onto principal component space\n\n");
	for(cluster = cluster_0, j = 0; j < cluster_0->numclusters; cluster = cluster->next, j++)
	{
		fprintf(fo, "cluster %d: water %d\n", j, cluster->wtrnum);
		fprintf(fo, "eigenvalue:\n");
		fprintf(fo, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", cluster->eigenvl[0], cluster->eigenvl[1], cluster->eigenvl[2], cluster->eigenvl[3], cluster->eigenvl[4], cluster->eigenvl[5]);
		fprintf(fo, "projected coor:\n");
		for(m = 0; m < cluster->wtrnum; m++)
			fprintf(fo, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", cluster->PCAcoor[0][m], cluster->PCAcoor[1][m], cluster->PCAcoor[2][m], cluster->PCAcoor[3][m], cluster->PCAcoor[4][m], cluster->PCAcoor[5][m]);
		fprintf(fo, "\n");
	}
	fclose(fo);

	ProbabilityDistributionFunction(wtr, cluster_0);
	SimpsonsIntegral(wtr, cluster_0, engy);
}

//====================================================================================================
// CovarianceMatrix: generate zero-mean covariance matrix for scaled positions and Euler angles (3x3)
//====================================================================================================
void CovarianceMatrix3x3(water *wtr, hcluster *cluster_0, energy *engy)
{
	hcluster	*cluster;
	float		**coor, **cov, coor_mean[3], **newcoor, newcov[3][3];
	float		**Eula, **Eulacov, Eula_mean[3], **newEula, newEulacov[3][3];
	int			m, n, s, i, j, frm, wtn, cl;
	float		egvl[3], e[3];
	FILE		*fo;

	printf("3x3 covariance matrix...\n");
	//go over each cluster 
	cov = (float **)calloc(3, sizeof(float*));
	Eulacov = (float **)calloc(3, sizeof(float*));
	for(i = 0; i < 3; i++)
	{
		cov[i] = (float*)calloc(3, sizeof(float));
		Eulacov[i] = (float*)calloc(3, sizeof(float));
	}
	for(cluster = cluster_0, cl = 0; cl < cluster_0->numclusters; cluster = cluster->next, cl++)
	{
		//assign coordinate and euler angle to coor
		coor = (float**)calloc(cluster->wtrnum, sizeof(float*));
		Eula = (float**)calloc(cluster->wtrnum, sizeof(float*));
		for(i = 0; i < cluster->wtrnum; i++)
		{
			coor[i] = (float*)calloc(3, sizeof(float));
			Eula[i] = (float*)calloc(3, sizeof(float));
		}
		for(i = 0; i < 3; i++)
		{
			coor_mean[i] = 0.0;
			Eula_mean[i] = 0.0;
		}
		for(m = 0; m < 3; m++)
		{
			for(n = 0; n < 3; n++)
			{
				cov[m][n] = 0.0;
				newcov[m][n]=0.0;
				Eulacov[m][n] = 0.0;
				newEulacov[m][n]=0.0;
			}
		}

		m = 0;
		for(frm = 0; frm < wtr->numfr; frm++)
		{
			if(cluster->wtrindex[frm] >= 0)
			{
				wtn = cluster->wtrindex[frm];
				for(i = 0; i < 3; i++)
				{
					coor[m][i] = wtr->O[frm][wtn].x[i];
					Eula[m][i] = cluster->wtrEulerAgl[frm][i];
				}
				m++;
			}
		}
		if(par.grid_energy == 0) printf("water inside cluster%d: %d (=%d)\n", cl, m, cluster->wtrnum);
		
		//calculate mean 
		for(m = 0; m < cluster->wtrnum; m++)
		{	
			for(i = 0; i < 3; i++)
			{
				coor_mean[i] += coor[m][i];
				Eula_mean[i] += Eula[m][i];
			}
		}		

		for(i = 0; i < 3; i++)
		{
			coor_mean[i] /= cluster->wtrnum;
			Eula_mean[i] /= cluster->wtrnum;
		}

		for(i = 0; i < cluster->wtrnum; i++)
		{
			for(m = 0; m < 3; m++)
			{
				coor[i][m] -= coor_mean[m];
				Eula[i][m] -= Eula_mean[m];
			}				
		}

		//construct covariance matrix
		for(m = 0; m < 3; m++)
		{
			for(n = 0; n < 3; n++)
			{
				for(i = 0; i < cluster->wtrnum; i++)
				{
					cov[m][n] += coor[i][m]*coor[i][n];
					Eulacov[m][n] += Eula[i][m]*Eula[i][n];
				}				
			}
		}

		for(m = 0; m < 3; m++)
		{
			for(n = 0; n < 3; n++)
			{
				if(cluster->wtrnum == 1)	
				{
					cov[m][n] /= cluster->wtrnum;
					Eulacov[m][n] /= cluster->wtrnum;
				}
				else						
				{
					cov[m][n] /= (cluster->wtrnum - 1);  // ???? n or n-1???
				//cluster->cov[m][n] = cov[m][n];
					Eulacov[m][n] /= (cluster->wtrnum - 1);  // ???? n or n-1???
				//cluster->Eulacov[m][n] = cov[m][n];
				}
			}
		}	
/*
		printf("covariance matrix A:\n");
		for(m = 0; m < 3; m++)
			printf("%10.6f %10.6f %10.6f\n", cov[m][0], cov[m][1], cov[m][2]);
		printf("covariance matrix B:\n");
		for(m = 0; m < 3; m++)
			printf("%10.6f %10.6f %10.6f\n", Eulacov[m][0], Eulacov[m][1], Eulacov[m][2]);
*/

		//calculate eigenvector and eigenvalue
		tred2(cov, 3, egvl, e);
		tqli(egvl, e, 3, cov);

/*		printf("eigenvectors A:\n");
		for(m = 0; m < 3; m++)
			printf("%10.6f %10.6f %10.6f: %f\n", cov[0][m], cov[1][m], cov[2][m], egvl[m]);
*/
		for(m = 0; m < 3; m++)
			cluster->eigenvl[m] = egvl[m];

		tred2(Eulacov, 3, egvl, e);
		tqli(egvl, e, 3, Eulacov);

/*		printf("eigenvectors B:\n");
		for(m = 0; m < 3; m++)
			printf("%10.6f %10.6f %10.6f: %f\n", Eulacov[0][m], Eulacov[1][m], Eulacov[2][m], egvl[m]);*/

		for(m = 0; m < 3; m++)
			cluster->eigenvl[m+3] = egvl[m];

		newcoor = (float**)calloc(3, sizeof(float*));
		newEula = (float**)calloc(3, sizeof(float*));
		for(m = 0; m < 3; m++)
		{
			newcoor[m] = (float*)calloc(cluster->wtrnum, sizeof(float));
			newEula[m] = (float*)calloc(cluster->wtrnum, sizeof(float));
			for(n = 0; n < cluster->wtrnum; n++)
			{
				newcoor[m][n] = 0.0;
				newEula[m][n] = 0.0;
			}
		}
		
		for(s = 0; s < cluster->wtrnum; s++)
		{
			for(m = 0; m < 3; m++)
			{
				for(n = 0; n < 3; n++)
				{
					newcoor[m][s] += cov[n][m] * coor[s][n];
					if(n == 2)
					{
						newEula[m][s] += Eulacov[n][m] * Eula[s][n];// * sin(Eula[s][n]);
					}
					else
					{
						newEula[m][s] += Eulacov[n][m] * Eula[s][n];
					}
				}
			}
		}

		cluster->PCAcoor = (float**)calloc(6, sizeof(float*));
		for(m = 0; m < 3; m++)
		{
			cluster->PCAcoor[m] = (float*)calloc(cluster->wtrnum, sizeof(float));
			cluster->PCAcoor[m+3] = (float*)calloc(cluster->wtrnum, sizeof(float));
			for(s = 0; s < cluster->wtrnum; s++)
			{
				cluster->PCAcoor[m][s] = newcoor[m][s];
				cluster->PCAcoor[m+3][s] = newEula[m][s];
			}
		}

		//free memory
		for(i = 0; i < cluster->wtrnum; i++)
		{
			free(coor[i]);
			free(Eula[i]);
		}
		free(coor);	
		free(Eula);
		for(i = 0; i < 3; i++)
		{
			free(newcoor[i]);
			free(newEula[i]);
		}
		free(newcoor);
		free(newEula);
	}
	
	for(m = 0; m < 3; m++)
	{
		free(cov[m]);
		free(Eulacov[m]);
	}
	free(cov);
	free(Eulacov);
	
	//-----------------OUTPUT file for histogram---------------//
	fo = fopen("PCAhistogram3x3.txt", "w");
	fprintf(fo, "This file has the coordinates projected onto principal component space\n\n");
	for(cluster = cluster_0, j = 0; j < cluster_0->numclusters; cluster = cluster->next, j++)
	{
		fprintf(fo, "cluster %d: water %d\n", j, cluster->wtrnum);
		fprintf(fo, "eigenvalue:\n");
		fprintf(fo, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", cluster->eigenvl[0], cluster->eigenvl[1], cluster->eigenvl[2], cluster->eigenvl[3], cluster->eigenvl[4], cluster->eigenvl[5]);
		fprintf(fo, "projected coor:\n");
		for(m = 0; m < cluster->wtrnum; m++)
			fprintf(fo, "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", cluster->PCAcoor[0][m], cluster->PCAcoor[1][m], cluster->PCAcoor[2][m], cluster->PCAcoor[3][m], cluster->PCAcoor[4][m], cluster->PCAcoor[5][m]);
		fprintf(fo, "\n");
	}
	fclose(fo);

	ProbabilityDistributionFunction(wtr, cluster_0);
	SimpsonsIntegral(wtr, cluster_0, engy);
}

//====================================================================================================
// ProbabilityDistributionFunction: calculate probability distribution function on each PC dimension
//====================================================================================================
void ProbabilityDistributionFunction(water *wtr, hcluster *cluster_0)
{ 
	int			i, j, cl, bin;
	hcluster	*cluster;
	float		n, min, max, interval, probability[par.numbins], fraction;
	float		integral, p;
	FILE		*fo;
	fo = fopen("PDF", "w");
	fprintf(fo, "Probability Density Function:\n\n");
	// go over each cluster 
	for(cluster = cluster_0, cl = 0; cl < cluster_0->numclusters; cluster = cluster->next, cl++)
	{
		fprintf(fo, "cluster: %d\n", cl);
		// go over each principle component dimension 
		for(i = 0; i < 6; i++)
		{
			fprintf(fo, "PC %d:", i);
			min = 9999.0;
			max = -9999.0;
			// find the min and max of this dimension 
			for(j = 0; j < cluster->wtrnum; j++)
			{
				if(cluster->PCAcoor[i][j] < min)
				{
					min = cluster->PCAcoor[i][j];
				}
				if(cluster->PCAcoor[i][j] > max)
				{
					max = cluster->PCAcoor[i][j];
				}
			}
			//max += 0.0001;
			//min -= 0.0001;
			interval = (max - min)/par.numbins;		// bin = 70
			fprintf(fo, "min max interval: %f %f %f\n", min, max, interval);
			cluster->min[i] = min;
			cluster->max[i] = max;
			cluster->interval[i] = interval;
			for(j = 0; j < par.numbins; j++)
			{
				probability[j] = 0.0;
			}
			
			for(j = 0; j < cluster->wtrnum; j++)
			{
				fraction = modff((cluster->PCAcoor[i][j] - min)/interval, &n);
				bin = (int)(n+Null);
				//printf("n :%d\n", bin);
				probability[bin]++;
			}
			
			fprintf(fo, "probability:\n");
			for(j = 0; j < par.numbins; j++)
			{
				probability[j] /= cluster->wtrnum*interval;
				cluster->probability[i][j] = probability[j];
				fprintf(fo, "%8f ", probability[j]);
				if((j+1)%10==0)	fprintf(fo, "\n");
			}
			fprintf(fo, "\n");
			
		}		
	}
	fclose(fo);	
}
//====================================================================================================
// SimpsonsIntegral: calculate probability density function on each PC dimension
//====================================================================================================
void SimpsonsIntegral(water *wtr, hcluster *cluster_0, energy *engy)
{
	int			i, j, cl;
	float		integral[6], p, lnp, entropy;
	hcluster	*cluster;
	FILE		*fo;

	engy->entropy = (float*) calloc(cluster_0->numclusters, sizeof(float));
	engy->grid_num = (int*) calloc(cluster_0->numclusters, sizeof(int));
	engy->coor = (float**) calloc(cluster_0->numclusters, sizeof(float*));
	engy->radius = (float*) calloc(cluster_0->numclusters, sizeof(float));
	engy->overall_occupancy = (float*)calloc(cluster_0->numclusters,sizeof(float));
	engy->average_occupancy = (float*)calloc(cluster_0->numclusters,sizeof(float));

	fo = fopen("entropy.egy", "w");
	fprintf(fo, "@cluster# entropy: entropy*Temp(298.15K) on each dimension (translation, rotation. Unit: kcal/mol)\n%4d\n", cluster_0->numclusters);
	// go over each cluster 
	for(cluster = cluster_0, cl = 0; cl < cluster_0->numclusters; cluster = cluster->next, cl++)
	{
		// go over each principle component dimension 
		entropy = 0.0;
		for(i = 0; i < 6; i++)
		{
			integral[i] = 0.0;
			
			for(j = 0; j < par.numbins; j++)
			{
				p = cluster->probability[i][j];
				lnp = log(p);
				if(p < Null)	lnp = 0.0;
				
				if(j==0 || j==par.numbins-1)
				{
					integral[i] += lnp*p;
				}
				else if(j%2 == 1)
				{
					integral[i] += 4*lnp*p;
				}
				else
				{
					integral[i] += 2*lnp*p;
				}
			}
			
			integral[i] *= -R*T0*cluster->interval[i]*0.001/3.0;		// unit: kcal/mol//S = -R*integral(p*ln(p))
			cluster->integral[i] = integral[i];
			entropy += integral[i];
		}
		//entropy += par.NeatS;	//refer to Minh2005 for the calculation of the neat water entropy

		entropy += BulkS;
		cluster->radius = 1.0;	//need to be calculated based on water occupancy and its distribution
		fprintf(fo, "%4d %10.3f : %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f    %7.4f %7.4f %7.4f    %7.2f    %7.4f %7.4f    %4d\n", cl, entropy, integral[0], integral[1], integral[2], integral[3], integral[4], integral[5], cluster->center[0], cluster->center[1], cluster->center[2], cluster->radius, cluster->overall_occupancy, cluster->average_occupancy, cluster->wtrnum);

        // YY 12.6.2016
        // Directly put TS as entropy
		//printf("%4d %10.3f : %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f    %7.4f %7.4f %7.4f    %7.2f    %7.4f %7.4f    %4d\n", cl, entropy, integral[0], integral[1], integral[2], integral[3], integral[4], integral[5], cluster->center[0], cluster->center[1], cluster->center[2], cluster->radius, cluster->overall_occupancy, cluster->average_occupancy, cluster->wtrnum);
        engy->grid_num[cl] = cluster->grid_num;                            //mind the sign!!!
		//printf("grid num: %d %d\n", cluster->grid_num, engy->grid_num[cl]);
        engy->entropy[cl] = entropy;                            //mind the sign!!!
        //engy->entropy[cl] = entropy*1000/T0;                  //mind the sign!!!
		engy->overall_occupancy[cl] = cluster->overall_occupancy;
		engy->average_occupancy[cl] = cluster->average_occupancy;
		engy->coor[cl] = (float*) calloc(3, sizeof(float));
        engy->coor[cl][0] = cluster->center[0];
        engy->coor[cl][1] = cluster->center[1];
        engy->coor[cl][2] = cluster->center[2];
        engy->radius[cl]  = cluster->radius;

	}	// YY entropy (TS)
	fclose(fo);
}
