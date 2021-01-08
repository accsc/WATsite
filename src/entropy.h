#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <sys/stat.h>
#include <unistd.h>

#include "molecule.h"
#include "parameters.h"
#include "util.h"

#ifndef _ENTROPY_H
#define _ENTROPY_H

void CalcEntropy(water *, hcluster *, energy *);

void EulerAngle(water *, hcluster *);
void CovarianceMatrix6x6(water *, hcluster *, energy *);
void CovarianceMatrix3x3(water *, hcluster *, energy *);
void ProbabilityDistributionFunction(water *, hcluster *);
void SimpsonsIntegral(water *, hcluster *, energy *);

#endif
