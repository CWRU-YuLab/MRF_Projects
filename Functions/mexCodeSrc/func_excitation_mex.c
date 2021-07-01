#include <mex.h> 
#include <stdio.h>
#include <math.h>

double calcrotmat(double nx, double ny, double *nz, int nPool, double *rotMtx)
{
    double ar, ai, br, bi, hp, cp, sp;
    double arar, aiai, arai2, brbr, bibi, brbi2, arbi2, aibr2, arbr2, aibi2;
    double phi;
    int cnt;
    int ii;
    
    for(ii=0; ii<(nPool*3*nPool*3); ii++)
        rotMtx[ii] = 0;
    
    for(ii=0; ii<nPool; ii++) 
    {
        cnt = nPool*3*3*ii+3*ii;
        if(nx==0 && ny==0 && nz[ii]==0)
        {
            rotMtx[0+cnt] = 1;
            rotMtx[nPool*3+1+cnt] = 1;
            rotMtx[nPool*3*2+2+cnt] = 1;
        }
        else
        {
            phi = sqrt(nx*nx+ny*ny+nz[ii]*nz[ii]);
            hp = phi/2;
            cp = cos(hp);
            sp = sin(hp)/phi;
            ar = cp;
            ai = -nz[ii]*sp;
            br = ny*sp;
            bi = -nx*sp;
            
            arar = ar*ar;
            aiai = ai*ai;
            arai2 = 2*ar*ai;
            brbr = br*br;
            bibi = bi*bi;
            brbi2 = 2*br*bi;
            arbi2 = 2*ar*bi;
            aibr2 = 2*ai*br;
            arbr2 = 2*ar*br;
            aibi2 = 2*ai*bi;
            
            rotMtx[0+cnt] = arar-aiai-brbr+bibi;
            rotMtx[1+cnt] = -arai2-brbi2;
            rotMtx[2+cnt] = -arbr2+aibi2;
            rotMtx[nPool*3+cnt] = arai2-brbi2;
            rotMtx[nPool*3+1+cnt] = arar-aiai+brbr-bibi;
            rotMtx[nPool*3+2+cnt] = -aibr2-arbi2;
            rotMtx[nPool*3*2+cnt] = arbr2+aibi2;
            rotMtx[nPool*3*2+1+cnt] = arbi2-aibr2;
            rotMtx[nPool*3*2+2+cnt] = arar+aiai-brbr-bibi;
        }
    }
}

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
    double *b1rDt;
    double *b1iDt;
    double *dfDt;
    double *decMtx;
    double *decVec;
    double *M0;
    int nPool;
    double *Mout;
    int ii, jj, nt;
        
    b1rDt = mxGetPr(prhs[0]);
    b1iDt = mxGetPr(prhs[1]);
    dfDt = mxGetPr(prhs[2]);
    decMtx = mxGetPr(prhs[3]);
    decVec = mxGetPr(prhs[4]);
    M0 = mxGetPr(prhs[5]);
    
    nPool = mxGetM(prhs[2]);
    nt = mxGetN(prhs[0]);
        
    double rotMtx[nPool*3*nPool*3];
    double vectemp1[nPool*3];
    double vectemp2[nPool*3];
    
    plhs[0] = mxCreateDoubleMatrix(nPool*3,1,mxREAL);
    Mout = mxGetPr(plhs[0]);
    
    for(ii=0; ii<nPool*3; ii++)
        vectemp1[ii] = *M0++;        
    
    for(ii=0; ii<nt; ii++)
    {
        calcrotmat(b1rDt[ii],b1iDt[ii],dfDt,nPool,rotMtx);
        multmatvec(rotMtx,vectemp1,vectemp2,nPool);
        multmatvec(decMtx,vectemp2,vectemp1,nPool);
        addvecs(vectemp1,decVec,vectemp2,nPool);
        for(jj=0;jj<nPool*3;jj++)
        {
            vectemp1[jj] = vectemp2[jj];
        }
    }
    for(ii=0;ii<nPool*3;ii++)
    {
        *Mout++ = vectemp1[ii];
    }
}