#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef _UTIL_H
#define _UTIL_H


static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
int *VectorInitialize(int, int);
float Distance_Point_Point(float *, float *);
void ResetVector(int *, int , int );
int Compare_Vectors(float *, float *, float);
void CreateNormVector(float *, float *, float *);
void NormalizeVector(float *);
void DistSq(float *, float *, float *);
void MakeVector(float *, float *, float *);
void VxV(float *,float *,float *);
void VdotV(float *,float *,float *);
float* MatrixMult(float [][3], float *);
void tred2(float **, int , float [], float []);
void tqli(float [], float [], int , float **);
float pythag(float , float);
void GetChars ( char [],int ,int ,char []);


#endif
