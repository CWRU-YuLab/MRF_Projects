#include <mex.h>
#include <stdio.h>
#include <math.h>

void multmatvec(double *mat, double *vec, double *matvec, int nPool)
{
    int i;
    int j;
    
    for (i=0; i<nPool*3; i++)
    {
        *matvec = 0.0;
        for (j=0; j<nPool*3; j++)
        {
            *matvec += mat[j*nPool*3+i]*vec[j];
        }
        matvec++;
    }
}

void addvecs(double *vec1, double *vec2, double *vecsum, int nPool)
{
    int i;
    for (i=0; i<nPool*3; i++)
        *vecsum++ = *vec1++ + *vec2++;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *rotMtx;
    double *decMtx;
    double *decVec;
    double *M0;
    double *Mout;
    int Nt;
    int nPool;
    int i;
    int j;
    
    decMtx = mxGetPr(prhs[0]);
    decVec = mxGetPr(prhs[1]);
    Nt = (int)(*mxGetPr(prhs[2]));
    M0 = mxGetPr(prhs[3]);
    nPool = mxGetM(prhs[3])/3;
    
    double vectemp1[nPool*3];
    double vectemp2[nPool*3];
    
    plhs[0] = mxCreateDoubleMatrix(nPool*3,Nt,mxREAL);
    Mout = mxGetPr(plhs[0]);
    
    for(i=0; i<nPool*3; i++)
    {
        vectemp1[i] = *M0++;
    }
    
    for(i=0; i<Nt; i++)
    {
        multmatvec(decMtx,vectemp1,vectemp2,nPool);
        addvecs(vectemp2,decVec,vectemp1,nPool);
        for(j=0; j<nPool*3; j++)
        {
            *Mout++ = vectemp1[j];
        }        
    }
}