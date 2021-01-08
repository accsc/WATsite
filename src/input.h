#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include "parameters.h"
#include "molecule.h"
#include "util.h"
#ifndef _INPUT_H
#define _INPUT_H

void ReadUsersParameters(int, char **);
void ReadCommandFile(char *);
void ReadBindingSite(char *, ligand *);
void DefineBindingPocket(ligand *, grid *);
void ReadEnthalpy(water *);
void ReadEnthalpy_fromCluster(water *);
void CutFrames();
void CleanCutFrames();

#endif
