#include "quadFunctions.h"
#include <cmath>

using namespace Ipopt;

const Number gearth = 9.81;
const Number m0 = 1.;

Number computeDynamicX(Number x1,Number vx,Number vx1,Number vxm,Number coeff){
    return x1 - coeff*(vx + vx1 + 4*vxm);
}

Number computeDynamicVX(Number vx1, Number u1, Number theta, Number u1_1, Number theta_1, Number u1m, Number thetam, Number coeff){
    return vx1 - coeff*(u1*sin(theta)/m0 + u1_1*sin(theta_1)/m0 + 4*(u1m*sin(thetam)/m0));
}

Number computeDynamicZ(Number z1,Number vz,Number vz1,Number vzm,Number coeff){
    return z1 - coeff*(vz + vz1 + 4 * vzm);
}

Number computeDynamicVZ(Number vz1,Number u1,Number theta,Number u1_1,Number theta1,Number u1m, Number thetam, Number coeff){
    return vz1 - coeff*((u1*cos(theta)/m0-gearth) + (u1_1*cos(theta1)/m0-gearth) + 4*(u1m*cos(thetam)/m0-gearth));
}

Number computeTheta(Number theta1,Number u2,Number u2_1,Number u2m, Number coeff){
    return theta1 - coeff*(u2+u2_1+4*u2m);
}

Number computeXM(Number x,Number x1,Number vx,Number vx1,Number coeff){    
    return (x+x1)/2 + coeff*(vx - vx1);
}

Number computeVXM(Number vx,Number vx1,Number u1,Number theta,Number u1_1,Number theta1,Number coeff){
    return (vx + vx1)/2 + coeff*(u1*sin(theta)/m0 - u1_1*sin(theta1)/m0);
}

Number computeZM(Number z,Number z1,Number vz,Number vz1,Number coeff){
    return (z + z1)/2 + coeff*(vz - vz1);
}

Number computeVZM(Number vz,Number vz1,Number u1,Number theta,Number u1_1,Number theta1, Number coeff){
    return (vz+vz1)/2 + coeff*((u1*cos(theta)/m0 - gearth) - (u1_1*cos(theta1)/m0 - gearth));
}

Number computeThetaM(Number theta,Number theta1,Number u2,Number u2_1, Number coeff){
    return (theta+theta1)/2 + coeff*(u2 - u2_1);
}

void getStateVectorConst(StateVectorConst * stVec,const Number * x, Index N){

    stVec->tf = x[N*7 + 2*(N-1)];

    for (int i = 0 ; i<N ; i++){
        stVec->x[i]      = x[i];
        stVec->vx[i]     = x[N+i];
        stVec->z[i]      = x[2*N+i];
        stVec->vz[i]     = x[3*N+i];
        stVec->theta[i]  = x[4*N+i];
        stVec->u1[i]     = x[5*N+i];
        stVec->u2[i]     = x[6*N+i];
    }

    Number coeff = stVec->tf/(N-1)/8;

    for (Index j = 0; j<N-1; j++){
        stVec->u1m[j] = x[7*N+j];
        stVec->u2m[j] = x[8*N-1+j];
        stVec->xm[j]    = computeXM(    stVec->x[j],    stVec->x[j+1],      stVec->vx[j], stVec->vx[j+1], coeff);
        stVec->vxm[j]   = computeVXM(   stVec->vx[j],   stVec->vx[j+1],     stVec->u1[j], stVec->theta[j],stVec->u1[j+1],stVec->theta[j+1], coeff);
        stVec->zm[j]    = computeZM(    stVec->z[j],    stVec->z[j+1],      stVec->vz[j], stVec->vz[j+1], coeff);
        stVec->vzm[j]   = computeVZM(   stVec->vz[j],   stVec->vz[j+1],     stVec->u1[j], stVec->theta[j],stVec->u1[j+1], stVec->theta[j+1], coeff);
        stVec->thetam[j] =computeThetaM(stVec->theta[j], stVec->theta[j+1], stVec->u2[j], stVec->u2[j+1], coeff);
    }
}

void getStateVector(StateVector * stVec, Number * x, Index N){

    stVec->tf = &x[9*N-2];

    for (Index i = 0 ; i<N ; i++){
        stVec->x[i] =   &x[i];
        stVec->vx[i] =  &x[N+i];
        stVec->z[i] =   &x[2*N+i];
        stVec->vz[i] =  &x[3*N+i];
        stVec->theta[i] = &x[4*N+i];
        stVec->u1[i] =  &x[5*N+i];
        stVec->u2[i] =  &x[6*N+i];
    }

    Number coeff = *stVec->tf/(N-1)/8;

    for (Index j = 0; j< N-1; j++){
        stVec->u1m[j] = &x[7*N+j];
        stVec->u2m[j] = &x[8*N-1+j];
        stVec->xm[j] = computeXM(   *stVec->x[j],       *stVec->x[j+1],         *stVec->vx[j], *stVec->vx[j+1],     coeff);
        stVec->vxm[j] = computeVXM( *stVec->vx[j],      *stVec->vx[j+1],        *stVec->u1[j], *stVec->theta[j], *stVec->u1[j+1], *stVec->theta[j+1], coeff);
        stVec->zm[j] = computeZM(   *stVec->z[j],       *stVec->z[j+1],         *stVec->vz[j], *stVec->vz[j+1],     coeff);
        stVec->vzm[j] = computeVZM( *stVec->vz[j],      *stVec->vz[j+1],        *stVec->u1[j], *stVec->theta[j], *stVec->u1[j+1], *stVec->theta[j+1], coeff);
        stVec->thetam[j] =computeThetaM(*stVec->theta[j], *stVec->theta[j+1],   *stVec->u2[j], *stVec->u2[j+1],     coeff);
    }
}

void getConstraintVector(ConstraintVector * cV, Number * g, Index N){

    for(int i = 0; i<5 ; i++){
        cV->boundCond[i] = &g[i];
        cV->finalCond[i] = &g[i+5];
    }

    Index upperB = N-1;
    Index ind = 0;

    for(int j = 10; j<upperB+10; j++){
        cV->dynamicX[ind] = &g[j];
        cV->dynamicVX[ind] = &g[j+upperB];
        cV->dynamicZ[ind] = &g[j+2*upperB];
        cV->dynamicVZ[ind] = &g[j+3*upperB];
        cV->dynamicTheta[ind] = &g[j+4*upperB];
        ind++;
    }
}

