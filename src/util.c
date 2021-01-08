#include "util.h"

/**********************************************************
    VectorInitialize
***********************************************************/
int *VectorInitialize(int k, int value)
{
  int   i;
  int   *t;
  if((t = (int *) calloc(k, sizeof(int))) == NULL)
  {
    printf("Insufficient Memory for allocating t.\n");
    exit( EXIT_FAILURE);
  }

  for(i = 0; i < k; i++)
  {
    t[i] = value;
  }
  return t;
}

/*********************************************************************************************
    Compare_Vectors
 *********************************************************************************************/
float Distance_Point_Point(float *a, float *b)
{
    return( sqrt(pow(a[0]-b[0], 2)+pow(a[1]-b[1], 2)+pow(a[2]-b[2], 2)) );
}

/**********************************************************
    ResetVector
***********************************************************/
void ResetVector(int *t, int k, int value)
{
  int i;

  for (i = 0; i < k; i++)
  {
    t[i] = value;
  }
}
/*********************************************************************************************
    Compare_Vectors
 *********************************************************************************************/
int Compare_Vectors(float *a, float *b, float value)
{
    if( fabs(a[0] - b[0]) < value )
    {
        if( fabs(a[1] - b[1]) < value )
        {
            if( fabs(a[2] - b[2]) < value )
            {
                return 0;
            }
            else if( a[2] < b[2] )
            {
                return -1;
            }
            else return 1;
        }
        else if (a[1] < b[1])
        {
            return -1;
        }
        else return 1;
    }
    else if (a[0] < b[0])
    {
        return -1;
    }
    else return 1;
}
//=================================================================================================
// CreateNormVector
//=================================================================================================
void CreateNormVector(float v1[], float v2[], float v3[])
{
    int         i;
    float       norm;

    norm  =  0.0;
    for(i = 0; i < 3; i++)
    {
        v3[i] = v1[i] - v2[i];
        norm += v3[i]*v3[i];
    }
    norm = sqrt(norm);
    if(norm == 0.0)
        return;

    for(i = 0; i < 3; i++)
        v3[i] /= norm;
}

//=================================================================================================
// NormalizeVector
//=================================================================================================
void NormalizeVector(float v[])
{
    int         i;
    float       norm;

    norm = 0.0;
    for(i = 0; i < 3; i++)
        norm += v[i]*v[i];
    norm = sqrt(norm);
    if(norm == 0.0)
        return;

    for(i = 0; i < 3; i++)
        v[i] /= norm;
}
//=================================================================================================
// DistSq
//=================================================================================================
void DistSq(float v1[], float v2[], float *dd)
{
    float d;

    *dd = 0.0;
    d = *v1 - (*v2);
    *dd = d*d;

    d = *(v1+1) - (*(v2+1));
    *dd += d*d;

    d = *(v1+2) - (*(v2+2));
    *dd += d*d;
}

//=================================================================================================
// MakeVector, make a normalized vector v3 from v1 and v2
//=================================================================================================
void MakeVector(float v1[], float v2[], float v3[])
{
    int         k;
    float       norm;

    norm  =  0.0;
    for(k = 0; k < 3; k++)
    {
        *(v3+k) = *(v1+k) - (*(v2+k));
        norm += *(v3+k)*(*(v3+k));
    }
    norm = sqrt(norm);
    if(norm == 0.0)
        return;

    for(k = 0; k < 3; k++)
        *(v3+k) = *(v3+k)/norm;
}

//=================================================================================================
// VxV
//=================================================================================================
void VxV(float v1[],float v2[],float v3[])
{
    int         k;
    float   norm;

    *v3    =*(v1+1)*(*(v2+2))-(*(v1+2))*(*(v2+1));
    *(v3+1)=*(v1+2)*(*v2)    -(*v1)    *(*(v2+2));
    *(v3+2)=*v1    *(*(v2+1))-(*(v1+1))*(*v2);

    norm=0.0;
    for(k=0; k<3; k++)
        norm += *(v3+k)*(*(v3+k));
    norm=sqrt(norm);
    if(norm==0.0)
        return;

    for(k=0; k<3; k++)
        *(v3+k)=*(v3+k)/norm;
}

//=================================================================================================
// VdotV
//=================================================================================================
void VdotV(float v1[],float v2[],float *cosphi)
{
    int i;

    *cosphi=0.0;
    for(i=0; i<3; i++)
        *cosphi += *(v1+i)*(*(v2+i));
}

//=================================================================================================
// MatrixMult
//=================================================================================================
float* MatrixMult(float M[][3], float v[])
{
    int         i, j;
    static float        z[3];

    for (i = 0; i < 3; i++)
    {
        z[i] = 0.0;
        for (j = 0; j < 3; j++)
        {
            z[i] += M[i][j]*v[j];
        }
    }

    return z;
}

/*********************************************************************************************
        tred2
 *********************************************************************************************/
// Householder reduction of a real, symmetric matrix a[0..n-1][0..n-1]. On output, a is replaced
// by the orthogonal matrix Q effecting the transformation. d[0..n-1] returns the diagonal elements
// of the tridiagonal matrix, and e[0..n-1] the off-diagonal elements, with e[0]=0. Several
// statements, as noted in comments, can be omitted if only eigenvalues are to be found, in which
// case a contains no useful information on output. Otherwise they are to be included.

void tred2(float **a, int n, float d[], float e[])
{
        int         l, k, j, i;
        float         scale, hh, h, g, f; 

        for (i = n-1; i >= 1; i--)
        {
                l = i-1;
                h = scale = 0.0;
                if(l > 0)
                {
                        for(k = 0; k <= l; k++)
                        {
                                scale += fabs(a[i][k]);
                        }

                        if(scale == 0.0)                                         // Skip transformation.
                        {
                                e[i] = a[i][l];
                        }
                        
                        else
                        {
                                for(k = 0; k <= l; k++)
                                {
                                        a[i][k] /= scale;                         // Use scaled a's for transformation.
                                        h += a[i][k]*a[i][k];                // Form sigma in h.
                                }
                                f = a[i][l];
                                g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
                                
                                e[i] = scale*g;
                                h -= f*g;                                                 // Now h is equation (11.2.4).
                                a[i][l] = f-g;                                         // Store u in the ith row of a.
                                f = 0.0;
                                for(j = 0;j <= l; j++)
                                {
                                        /* Next statement can be omitted if eigenvectors not wanted */
                                        a[j][i] = a[i][j]/h;                 // Store u=H in ith column of a.
                                        g = 0.0;                                         // Form an element of A  u in g.
                                        for(k = 0; k <= j; k++)
                                        {
                                                g += a[j][k]*a[i][k];
                                        }
                                        for(k = j+1; k <= l; k++)
                                        {
                                                g += a[k][j]*a[i][k];
                                        }
                                        e[j] = g/h;                                 // Form element of p in temporarily unused element of e.
                                        f += e[j]*a[i][j];
                                }
                                hh = f/(h+h);                                         // Form K, equation (11.2.11).
                                for(j = 0;j <= l; j++)
                                {                                                                 // Form q and store in e overwriting p.
                                        f = a[i][j];
                                        e[j] = g = e[j]-hh*f;
                                        for(k = 0; k <= j; k++)         // Reduce a, equation (11.2.13).
                                        {
                                                a[j][k] -= (f*e[k]+g*a[i][k]);
                                        }
                                }
                        }
                }
                
                else
                {
                        e[i] = a[i][l];
                }
                d[i] = h;
        }
        /* Next statement can be omitted if eigenvectors not wanted */
        d[0] = 0.0;
        e[0] = 0.0;
        /* Contents of this loop can be omitted if eigenvectors not wanted except for statement d[i]=a[i][i]; */
        for(i = 0; i <= n-1; i++)
        {                                                                                         // Begin accumulation of transformationmatrices.
                l = i-1;
                if(d[i])
                {                                                                                 // This block skipped when i=1.
                        for(j = 0;j <= l; j++)
                        {
                                g = 0.0;
                                for (k = 0; k <= l; k++)                                 // Use u and u=H stored in a to form PQ.
                                {
                                        g += a[i][k]*a[k][j];
                                }
                                for(k = 0; k <= l; k++)
                                {
                                        a[k][j] -= g*a[k][i];
                                }
                        }
                }
                d[i] = a[i][i];                                                 // This statement remains.
                a[i][i] = 1.0;                                                         // Reset row and column of a to identity matrix for next iteration.
                for(j = 0; j <= l; j++)
                {
                        a[j][i] = a[i][j] = 0.0;
                }
        }
}



/*********************************************************************************************
        tqli
 *********************************************************************************************/
// QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of a real, symmetric,
// tridiagonal matrix, or of a real, symmetric matrix previously reduced by tred2 x11.2. On
// input, d[0..n-1] contains the diagonal elements of the tridiagonal matrix. On output, it returns
// the eigenvalues. The vector e[0..n-1] inputs the subdiagonal elements of the tridiagonal matrix,
// with e[0] arbitrary. On output e is destroyed. When finnding only the eigenvalues, several lines
// may be omitted, as noted in the comments. If the eigenvectors of a tridiagonal matrix are desired,
// the matrix z[0..n-1][0..n-1] is input as the identity matrix. If the eigenvectors of a matrix
// that has been reduced by tred2 are required, then z is input as the matrix output by tred2.
// In either case, the kth column of z returns the normalized eigenvector corresponding to d[k].

void tqli(float d[], float e[], int n, float **z)
{
//        float         pythag(float a, float b);
        int         m, l, iter, i, k;
        float         s, r, p, g, f, dd, c, b;

        for(i = 1; i <= n-1; i++)
        {
                e[i-1] = e[i];                                                                 // Convenient to renumber the elements of e.
        }
        e[n-1] = 0.0;
        for(l = 0; l <= n-1; l++)
        {
                iter = 0;
                do
                {
                        for(m = l; m <= n-2; m++)
                        {                                                                                 // Look for a single small subdiagonal element to split        the matrix.
                                dd = fabs(d[m])+fabs(d[m+1]);
                                if((float)(fabs(e[m])+dd) == dd) break;
                        }
                        if(m != l)
                        {
                                if(iter++ == 10000)
                                {
                                        //nrerror("Too many iterations in tqli");
                                        printf("Too many iterations in tqli\n");
                                }
                                g = (d[l+1]-d[l])/(2.0*e[l]);                 // Form shift.
                                r = pythag(g, 1.0);
                                g = d[m]-d[l]+e[l]/(g+SIGN(r,g));         // This is dm - ks.
                                s = c = 1.0;
                                p = 0.0;
                                for(i = m-1; i>=l; i--)
                                {                                                                        // A plane rotation as in the original QL, followed by Givens rotations to restore tridiagonal form.
                                        f = s*e[i];
                                        b = c*e[i];
                                        e[i+1] = (r = pythag(f,g));
                                        if(r == 0.0)
                                        {                                                                 // Recover from underflow.
                                                d[i+1] -= p;
                                                e[m] = 0.0;
                                                break;
                                        }
                                        s = f/r;
                                        c = g/r;
                                        g = d[i+1]-p;
                                        r = (d[i]-g)*s+2.0*c*b;
                                        d[i+1] = g+(p=s*r);
                                        g = c*r-b;
                                        /* Next loop can be omitted if eigenvectors not wanted*/
                                        for(k = 0; k <= n-1; k++)
                                        {                                                                 // Form eigenvectors.
                                                f = z[k][i+1];
                                                z[k][i+1] = s*z[k][i]+c*f;
                                                z[k][i] = c*z[k][i]-s*f;
                                        }
                                }
                                if(r == 0.0 && i >= l) continue;
                                d[l] -= p;
                                e[l] = g;
                                e[m] = 0.0;
                        }
                } while (m != l);
        }
}



/*********************************************************************************************
        pythag
 *********************************************************************************************/
// Computes (a2 + b2)^1/2 without destructive underflow or overflow.

float pythag(float a, float b)
{
        float         absa, absb;

        absa = fabs(a);
        absb = fabs(b);
        if(absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
        else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

/*********************************************************************************************
   GetChars  copies chars from n1 to n2 to res
*********************************************************************************************/
void GetChars ( char line[],int n1,int n2,char res[] )
{
   int    i,j;

   j=0;
   for ( i=n1; i<=n2; i++ )
       if ( line[i] != ' ' )
           res[j++]=line[i];
   res[j]='\0';
}



