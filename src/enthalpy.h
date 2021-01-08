
#ifndef _ENTHALPY_H
#define _ENTHALPY_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "parameters.h"
#include "molecule.h"

void CalculateFreeEnergy(water *, hcluster *, energy *);
void OutputClusterInfo(grid *, energy *);

#endif
