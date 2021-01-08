/***************************************************************************
 *       This script
 *       1) reads user parameters
 *   2) reads command file
 *   3) reads ligand file and defines the binding pocket                           *
 *       4) defines (grids) extending from the minimum and maximum point of ligand;                *
 ***************************************************************************/
#include "input.h"

//=================================================================================================
// ReadUsersParameters
//=================================================================================================
void ReadUsersParameters(int argc, char *argv[]) {
	char filecomm[kString];
	int j;
	// default:
	sprintf(filecomm, "input.bcf");
	sprintf(par.outfolder, "OUTPUT");

	// Read command line arguments
	// Standard name of command file

	for (j = 1; j < argc; j += 2) {
		if (!strcmp(argv[j], "-c")) {
			strcpy(filecomm, argv[j + 1]);
		} else if (!strcmp(argv[j], "-O") || !strcmp(argv[j], "-o")) {
			strcpy(par.outfolder, argv[j + 1]);
		} else if (!strcmp(argv[j], "-H") || !strcmp(argv[j], "-h")) {
			printf(
					"Usuage: watsite \n  -c       input.bcf \n  -o       <output folder> \n");
			exit(0);
		}
	}

	//read command file
	printf("command file %s\n", filecomm);
	ReadCommandFile(filecomm);
}

//=================================================================================================
// ReadCommandFile
//=================================================================================================
void ReadCommandFile(char *filename) {
	FILE *fp, *fi;
	char line[lString];

	fp = fopen(filename, "r");
	if (!fp) {
		printf("Command file %s is missing.\n", filename);
		exit(0);
	} else {
		// input of ligands
		while (!feof(fp)) {
			fgets(line, lString, fp);
			if (strstr(line, "$Input_files:"))
				break;
		}
		fgets(line, lString, fp);
		sscanf(line, "%s", par.LigandFile);
		printf("ligand file: %s\n", par.LigandFile);
		fgets(line, lString, fp);
		sscanf(line, "%s", par.ProtFile);
		printf("protein file: %s\n", par.ProtFile);
		fgets(line, lString, fp);
		sscanf(line, "%s", par.WaterEnthalpyFolder);
		printf("water enthalpy folder: %s\n", par.WaterEnthalpyFolder);

		rewind(fp);
		while (!feof(fp)) {
			fgets(line, lString, fp);
			if (strstr(line, "$Input_MDsnapshots:"))
				break;
		}
		fgets(line, lString, fp);
		sscanf(line, "%d", &par.numfr);
		printf("number of frames: %d\n", par.numfr);
		fgets(line, lString, fp);
		sscanf(line, "%s", par.water_model);
		printf("water model: %s\n", par.water_model);

		if (strcmp(par.water_model, "SPC/E") == 0) {
			par.site3_water = 1;
			par.E_bulk = BulkE_spce;
		} else if (strcmp(par.water_model, "TIP3P") == 0) {
			par.site3_water = 1;
			par.E_bulk = BulkE_tip3p;
		} else if (strcmp(par.water_model, "OPC") == 0) {
			par.site4_water = 1;
			par.E_bulk = BulkE_opc;
		} else if (strcmp(par.water_model, "TIP4P") == 0) {
			par.site4_water = 1;
			par.E_bulk = BulkE_tip4p;
		} else if (strcmp(par.water_model, "TIP4PEW") == 0) {
			par.site4_water = 1;
			par.E_bulk = BulkE_tip4pew;
		} else if (strcmp(par.water_model, "TIP5P") == 0) {
			par.site5_water = 1;
			par.E_bulk = BulkE_tip5p;
		} else if (strcmp(par.water_model, "AMOEBA") == 0) {
			par.site3_water = 1;
			par.E_bulk = BulkE_amoeba;
		} else {
			printf("Water Model %s is not supported...\n", par.water_model);
			exit(0);
		}

		rewind(fp);
		while (!feof(fp)) {
			fgets(line, kString, fp);
			if (strstr(line, "$Grid_Parameters:"))
				break;
		}
		fgets(line, lString, fp);
		sscanf(line, "%f", &(par.size_water));
		fgets(line, lString, fp);
		sscanf(line, "%f", &(par.griddelta));
		fgets(line, lString, fp);
		sscanf(line, "%f", &(par.FarestDist));
		fgets(line, lString, fp);
		sscanf(line, "%f", &(par.griddensity));
		printf("griddelta %f\n", par.griddelta);

		rewind(fp);
		while (!feof(fp)) {
			fgets(line, kString, fp);
			if (strstr(line, "$Cluster_method:"))
				break;
		}
		fgets(line, lString, fp);
		sscanf(line, "%d", &(par.cluster_method));
		fgets(line, lString, fp);
		sscanf(line, "%d", &(par.grid_energy));

		rewind(fp);
		while (!feof(fp)) {
			fgets(line, lString, fp);
			if (strstr(line, "$Hydro_Cluster:"))
				break;
		}
		fgets(line, lString, fp);
		sscanf(line, "%d", &(par.number_clusters));
		fgets(line, lString, fp);
		sscanf(line, "%f", &(par.cluster_mean_dist));
		fgets(line, lString, fp);
		sscanf(line, "%f", &(par.maxdist));
		fgets(line, lString, fp);
		sscanf(line, "%f", &(par.QT_occupancy_cutoff));

		rewind(fp);
		while (!feof(fp)) {
			fgets(line, kString, fp);
			if (strstr(line, "$DBSCAN:"))
				break;
		}
		fgets(line, lString, fp);
		sscanf(line, "%d", &(par.dbscan_START));
		fgets(line, lString, fp);
		sscanf(line, "%d", &(par.dbscan_END));

		rewind(fp);
		while (!feof(fp)) {
			fgets(line, kString, fp);
			if (strstr(line, "$Entropy:"))
				break;
		}
		fgets(line, lString, fp);
		sscanf(line, "%d", &(par.covdimension));
		fgets(line, lString, fp);
		sscanf(line, "%d", &(par.numbins));


		//MM 1/11/2018 New input parameters to use clustered trajectory
		rewind(fp);
		while(!feof(fp))
		{
			fgets(line, kString, fp);
			if(strstr(line, "$Clustered_Traj:")){
				par.clusteredTrajFlag = 1;
				break;
			}
		}

		if(par.clusteredTrajFlag == 1){
			fgets(line, lString, fp);
			sscanf(line, "%s", par.TrajEgyAssocFile);
			printf("traj-watEgy association file: %s\n", par.TrajEgyAssocFile);

			fgets(line, lString, fp);
			sscanf(line, "%d", &par.numclust);
			printf("cluster number: %d\n", par.numclust);

			fgets(line, lString, fp);
			sscanf(line, "%d", &par.totalnumfr);
			printf("total original frames: %d\n", par.totalnumfr);
		}

		/*
		//MM 1/11/2018 New input parameters to test convergence
		rewind(fp);
		while(!feof(fp))
		{
			fgets(line, kString, fp);
			if(strstr(line, "$Cut_Frame:")){
				par.frameCutoffFlag = 1;
				break;
			}
		}

		if(par.frameCutoffFlag == 1){
			fgets(line, lString, fp);
			sscanf(line, "%d", &par.frameCutoff);
			printf("frame cutoff: %d\n", par.frameCutoff);

			CutFrames();
		}
		*/
	}
	fclose(fp);
}

//=================================================================================================
// ReadUsersParameters
//=================================================================================================
void ReadBindingSite(char *infile, ligand *lig) {
	FILE *fi;
	int nl = 0, nh = 0, i, j;   // nh: number of heavy atoms
	char line[lString], atmn[8], coor[kString];
	float x, y, z;

	fi = fopen(infile, "r");
	if (!fi) {
		printf("ligand file %s is missing.\n", infile);
		exit(0);
	}

	else {
		rewind(fi);
		while (!feof(fi)) {
			fgets(line, lString, fi);
			if (strstr(line, "ATOM") || strstr(line, "HETATM")) {
				nl++;
				GetChars(line, 75, 81, atmn);
				if (!strstr(atmn, "H"))
					nh++;
			}

		}
		lig->atm = (atom *) calloc(nh, sizeof(atom));
		printf("number of heavy atoms in ligand: %d\n", nh);
		lig->numheavyatms = nh;

		rewind(fi);
		j = 0;
		while (!feof(fi)) {
			fgets(line, lString, fi);
			if (strstr(line, "ATOM") || strstr(line, "HETATM")) {
				GetChars(line, 75, 81, atmn);
				if (!strstr(atmn, "H")) {
					GetChars(line, 30, 37, coor);
					x = atof(coor);
					GetChars(line, 38, 45, coor);
					y = atof(coor);
					GetChars(line, 46, 54, coor);
					z = atof(coor);
					strcpy(lig->atm[j].name, atmn);
					lig->atm[j].x[0] = x;
					lig->atm[j].x[1] = y;
					lig->atm[j].x[2] = z;
					lig->atm[j].Num = j;
					j++;
					printf("atmn %sx y z: %f %f %f\n", atmn, x, y, z);
				}
			}
		}
	}
	printf("Finish reading ligand\n");
	fclose(fi);
}

//=================================================================================================
// DefineBindingPocket
//=================================================================================================
void DefineBindingPocket(ligand *lig, grid *grd) {
	int i, j, k, l, m, num;
	float longest, y[3];
	FILE *fo;

	for (j = 0; j < 3; j++) {
		grd->xmin[j] = 99999.9;
		grd->xmax[j] = -99999.9;
	}

	// determine size of binding pocket
	for (k = 0; k < lig->numheavyatms; k++) {
		for (j = 0; j < 3; j++) {
			if (lig->atm[k].x[j] < grd->xmin[j]) {
				grd->xmin[j] = lig->atm[k].x[j];
			}
			if (lig->atm[k].x[j] > grd->xmax[j]) {
				grd->xmax[j] = lig->atm[k].x[j];
			}
		}
	}

	// FarestDist is user-defined parameter (e.g. 10A)
	// griddelta is grid spacing (e.g. 0.25A)
	printf("binding site min %f %f %f\nbinding site max %f %f %f\n", grd->xmin[0], grd->xmin[1],
			grd->xmin[2], grd->xmax[0], grd->xmax[1], grd->xmax[2]);
	for (j = 0; j < 3; j++) {
		grd->xmin[j] -= par.FarestDist; // + 2; //YY 12/5/16 Only extend to user defined box size
		grd->xmax[j] += par.FarestDist;    // + 2;
		longest = grd->xmax[j] - grd->xmin[j];
		grd->GridSize[j] = par.griddelta * floor(longest / par.griddelta);
		printf("Axix %d: extending %f angstrom; longest %f\n", j, par.FarestDist, longest);
	}
	grd->griddelta = par.griddelta;
	// Determine zeropoint of Grid
	for (j = 0; j < 3; j++) {
		grd->ZeroPoint[j] = grd->xmin[j];
	}
	printf("ZeroPoint %f %f %f \nMaxPoint %f %f %f \nGridSize %f %f %f\n",
			grd->ZeroPoint[0], grd->ZeroPoint[1], grd->ZeroPoint[2],
			grd->xmax[0], grd->xmax[1], grd->xmax[2],
			grd->GridSize[0], grd->GridSize[1], grd->GridSize[2]);

	// Number of grid points
	for (j = 0; j < 3; j++) {
		grd->NumGP_r[j] = (int) (grd->GridSize[j] / par.griddelta) + 1;
	}
	printf("Number of grid points: %f %f %f %d %d %d\n",
			grd->GridSize[0] / par.griddelta, grd->GridSize[1] / par.griddelta,
			grd->GridSize[2] / par.griddelta, grd->NumGP_r[0], grd->NumGP_r[1], grd->NumGP_r[2]);
	grd->numGP1 = grd->NumGP_r[2]; //Z-axis
	grd->numGP2 = grd->NumGP_r[1] * grd->NumGP_r[2]; //Z-axis * Y-axis
	grd->numGP3 = grd->NumGP_r[0] * grd->NumGP_r[1] * grd->NumGP_r[2]; // 3 dimensions?

	// Coordinates of grid points
	grd->GP_r = (float **) calloc(grd->numGP3, sizeof(float *));
	for (k = 0; k < grd->NumGP_r[0]; k++) //YY NumGP_r[0] = the number of grid points on x axis
	{
		y[0] = k * grd->griddelta + grd->ZeroPoint[0];
		for (l = 0; l < grd->NumGP_r[1]; l++) {
			y[1] = l * grd->griddelta + grd->ZeroPoint[1];
			for (m = 0; m < grd->NumGP_r[2]; m++) {
				y[2] = m * grd->griddelta + grd->ZeroPoint[2];
				num = k * grd->numGP2 + l * grd->numGP1 + m; //YY num = index of the grid point

				grd->GP_r[num] = (float *) calloc(3, sizeof(float));

				grd->GP_r[num][0] = y[0];
				grd->GP_r[num][1] = y[1];
				grd->GP_r[num][2] = y[2];
				//YY store the coordinates of grid points: grd->GP_r[num][0-2] x,y,x coordinates of the num_th grid point
			}
		}
	}

	//YY_8-29-14: GridSize (grd->GridSize[j]) = distance; grd->GridSize[j]/par.griddelta = number of grid point;  grd->numGP3 = 3D number of grid points.

	printf("numGP %d = %d * %d * %d\n", grd->numGP3, grd->NumGP_r[0],
			grd->NumGP_r[1], grd->NumGP_r[2]);
	printf("Binding Pocket Defined!\n");

}

//=================================================================================================
// ReadWaterEnergies
//=================================================================================================
void ReadEnthalpy(water *wtr) {
	FILE *fi;
	char infile[kString];
	char line[lString];
	int step, interval, frame;
	int watnum;
	float watenergy;
	int numfr = 0, numatm = 0, numWAT = 0, watinside = 0;
	int dd, WATatm = 0;
	int i, j, k;

	numWAT = wtr->numWATs;
	numatm = wtr->numatms;
	//numfr  = wtr->numfr;

	int file_count = 0;
	DIR * dirin;
	struct dirent * entry;

	if ((dirin = opendir(par.WaterEnthalpyFolder)) == NULL) {
		fprintf(stderr, "Cannot Open Water Enthalpy Folder: %s\n",
				par.WaterEnthalpyFolder);
		exit(0);
	}
	while ((entry = readdir(dirin)) != NULL) {
		if ((!strcmp(entry->d_name, ".")) || (!strcmp(entry->d_name, ".."))) {
		} else {
			file_count++;
			//printf("%c  %s\n", entry->d_type, entry->d_name);
		}
	}
	closedir(dirin);
	//printf("Closed directory\n");

	if (file_count < wtr->numfr) {
		printf(
				"Error: Number of frames (%d) is more than the number of water enthalpy files (%d)\n",
				wtr->numfr, file_count);
		exit(0);
	} else {
		numfr = wtr->numfr;
	}
	//printf("Number of frames matches: %d %d \n", wtr->numfr, file_count);

	wtr->dH = (float**) calloc(numfr, sizeof(float*));
	//printf("Allocated dH\n");
	/*
	 go over water energy file for each frame of MD simulation
	 WaterEnergy/watEgy_1 ... watEgy_5000
	 store water energies into wat->dH[frame][waterindex in binding pocket]
	 */
	for (i = 0; i < numfr; i++) {
		//sprintf(infile, "%s/egy_%d",par.WaterEnthalpyFolder, i+1);
		sprintf(infile, "%s/watEgy_%d", par.WaterEnthalpyFolder, i + 1);
		fi = fopen(infile, "r");
		//printf("open water energy file: %s\n", infile);

		fgets(line, lString, fi);
		//printf("%s", line);

		if (strstr(line, "$STEP")) {
			sscanf(line, "%*s%d", &step);
		}

		watinside = wtr->WATinside[i] + 1;
		wtr->dH[i] = (float*) calloc(watinside, sizeof(float));

		for (j = 0; j < numWAT; j++) {
			//printf("frame i: %d  water index j: %d \n", i, j);
			fgets(line, lString, fi);
			sscanf(line, "%d%f", &watnum, &watenergy);
			//printf("Step: %d  Water Atom Number:  %d\n", step, watnum);

			if (watenergy > -0.000001 && watenergy < 0.000001) {
				printf(
						"Warning:\nsome of the water enthalpy values are incorrect!\nCheck the water energy file: %s",
						infile);
				printf("Step:%d  Water Number:%d Energy:%f\n", step, watnum, watenergy);
				exit(0);
			}
			for (k = 0; k < watinside; k++) {
				//printf("watnum in Egy (%d) -- atn in traj (%d)\n",watnum,wtr->atn[i][k]);
				if (watnum == wtr->atn[i][k]) {
					wtr->dH[i][k] = watenergy * 2;
					//printf("Frame(%d) Enthalpy: %f\n",i,watenergy);
				}
			}
		}
		fclose(fi);
	}
	printf("Finish reading water energies...\n");
}

//=================================================================================================
// ReadWaterEnergies from clustering
//=================================================================================================
void ReadEnthalpy_fromCluster(water *wtr)
{
	FILE    *fFrames, *fi;
	char    infile[kString], framefile[lString];
	char    line[lString];
	char 	buffer[10];
	int 	numclust;
	int     step, interval, frame, cluster;
	int     watnum;
	float   watenergy;
	int     numfr = 0, numatm = 0, numWAT = 0, watinside = 0, current_frame = 0, totalnumfr = 0;
	int     dd, WATatm=0;
	int     i, j, k;

	numWAT = wtr->numWATs;
	numatm = wtr->numatms;
	numfr  = wtr->numfr;
	totalnumfr = par.totalnumfr;

	int file_count = 0;
	DIR * dirin;
	struct dirent * entry;

	if ((dirin = opendir(par.WaterEnthalpyFolder)) == NULL) 
	{
		fprintf(stderr, "Cannot Open Water Enthalpy Folder: %s\n", par.WaterEnthalpyFolder);
		exit(0);
	}
	while ((entry = readdir(dirin)) != NULL) 
	{
		if (  (!strcmp(entry->d_name, ".") ) || (!strcmp(entry->d_name, ".."))  ) {
		}
		else
		{
			file_count++;
		}
	}
	closedir(dirin);

	if (file_count != totalnumfr)
	{
		printf("Error: Number of frames (%d) does not match with the number of water enthalpy files (%d)\n", totalnumfr, file_count);
		exit(0);
	}

	fFrames = fopen(par.TrajEgyAssocFile, "r");
	
	if(!fFrames)
	{
		printf("Cluster-frame association file %s is missing.\n", par.TrajEgyAssocFile);
		exit(0);
	}

	numclust = par.numclust;

	int frame_count = 0;
	rewind(fFrames);
	fgets(line, lString, fFrames); //MM move past header line
	while(!feof(fFrames))
	{
		fgets(line, lString, fFrames);
		GetChars(line, 11, 20, buffer);
		cluster = atoi(buffer);
		if(cluster == numclust)
		{
			frame_count++;
		}
	}

	int frames[frame_count];
	frame_count = 0;
	rewind(fFrames);
	fgets(line, lString, fFrames); //MM move past header line
	while(!feof(fFrames))
	{
		fgets(line, lString, fFrames);
		GetChars(line, 11, 20, buffer);
		cluster = atoi(buffer);
		if(cluster == numclust)
		{
			GetChars(line, 0, 7, buffer);
			frame = atoi(buffer);
			frames[frame_count] = frame;
			frame_count++;
			printf("frame count: %d\nframe: %d\n", frame_count, frame);
		}
	}
	if (par.frameCutoffFlag == 0 && frame_count != numfr){
		printf("Error: Number of frames (%d) in cluster (%d) does not match with the number of frames inputed (%d)\n", frame_count, numclust, numfr);
		exit(0);
	}

	wtr->dH = (float**)calloc(numfr, sizeof(float*));

	for(i = 0; i < numfr; i++)
	{
		current_frame = frames[i];
		sprintf(infile, "%s/watEgy_%d",par.WaterEnthalpyFolder, current_frame);	
		fi = fopen(infile, "r");
		printf("open water energy file: %s\n", infile);

		fgets(line, lString, fi);
		//printf("%s", line);
		if(strstr(line, "$STEP"))
		{
			sscanf( line, "%*s%d", &step );
		}
		watinside = wtr->WATinside[i] + 1;
		wtr->dH[i] = (float*)calloc(watinside, sizeof(float));

		for (j = 0; j < numWAT; j ++)
		{
			fgets(line, lString, fi);
			sscanf(line, "%d%f", &watnum, &watenergy);

			if ( watenergy > -0.000001 && watenergy < 0.000001 )
			{
                printf("Warning:\nsome of the water enthalpy values are incorrect!\nCheck the water energy file: %s", infile);
                printf("Step:%d  Water Number:%d Energy:%f\n", step, watnum, watenergy);
				//exit(0);
			}
			for (k = 0; k < watinside; k++)
			{
				if(watnum == wtr->atn[i][k])
				{
					wtr->dH[i][k] = watenergy*2;
					//printf("Frame(%d) Enthalpy: %f\n",i,watenergy);
				}
			}
		}
		fclose(fi);
	}
	printf("Finish reading water energies...\n");
}

//=================================================================================================
// Removes frames after n in trajectory			MM 1/18/18
//=================================================================================================
void CutFrames()
{
	FILE    *fi, *fo;
	char    line[lString];
	int firstFrames = par.frameCutoff;
	int n;
	n = 0;

	fi = fopen(par.ProtFile, "r");
	fo = fopen("/tmp/cuttraj.pdb", "w");
	if (!fi)
	{
		printf("Protein file %s is missing.\n", par.ProtFile);
		exit(0);
	}else{
		while(!feof(fi) && n < firstFrames){
			fgets(line, lString, fi);
			fputs(line, fo);
			if(strstr(line, "END")){
				n++;
			}
		}
		strcpy(par.ProtFile, "/tmp/cuttraj.pdb");
		par.numfr = par.frameCutoff;
	}
}

void CleanCutFrames()
{
	if (remove(par.ProtFile) == 0)
      printf("Temporary trajectory deleted successfully");
  	else
      printf("Unable to delete the temporary trajectory (%s)", par.ProtFile);
}