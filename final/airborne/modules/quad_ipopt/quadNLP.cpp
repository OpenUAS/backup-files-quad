#include "quadNLP.h"

#include <cassert>
#include <cstdio>
#include <cmath>

using namespace Ipopt;

// Set the constant parameters for the model
const Number maxthrust = 20.;
const Number minthrust = 1.;
const Number maxthetadot = 5.;
const Number m0 = 1.;
const Number gearth = 9.81;

const Number vxn = 0.;
const Number vzn = 0.;
const Number tn = 1.;

quadNLP::quadNLP(Index N, Number x0, Number z0, Number vx0, Number vz0, Number theta0, Number xn, Number zn, Number thetan, Number u10, Number u20) : N_(N) {
    this->x0 = x0;
    this->z0 = z0;
    this->vx0 = vx0;
    this->vz0 = vz0;
    this->theta0 = theta0;
    this->xn = xn;
    this->zn = zn;
    this->thetan = thetan;
    this->u10 = u10;
    this->u20 = u20;
}

quadNLP::~quadNLP() {}

// returns the size of the problem
bool quadNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                           Index& nnz_h_lag,
                           IndexStyleEnum& index_style)
{
    //Number of variables
    n = 9*N_-1;

    //Number of constraints
    m = 5*(N_-1) + 10;

    //Number of non-zeros in Jaccobian
    //    nnz_jac_g = m*n;
    nnz_jac_g = 44*(N_-1)+10;

    //Non-zeros in hessian(we use approximation
    nnz_h_lag = 0;

    // Use the C style indexing (0-based) for the matrices
    index_style = TNLP::C_STYLE;

    return true;
}

bool quadNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                              Index m, Number* g_l, Number* g_u)
{
    //State vector for lower and upper boundaries of variables
    StateVector * vecLow = new StateVector(N_);
    StateVector * vecUpp = new StateVector(N_);

    //Get addresses for variables
    getStateVector(vecLow,x_l,N_);
    getStateVector(vecUpp,x_u,N_);

    //Boundary for t final
    *vecLow->tf = 0.;
    *vecUpp->tf = 2e19;

    //Set boundaries for variables with dimension N
    for (int i = 0; i < N_; i++){
        *vecLow->x[i] = -2e19;
        *vecUpp->x[i] = 2e19;

        *vecLow->vx[i] = -2e19;
        *vecUpp->vx[i] = 2e19;

        *vecLow->z[i] = -2e19;
        *vecUpp->z[i] = 2e19;

        *vecLow->vz[i] = -2e19;
        *vecUpp->vz[i] = 2e19;

        *vecLow->theta[i] = -2e19;
        *vecUpp->theta[i] = 2e19;

        *vecLow->u1[i] = minthrust;
        *vecUpp->u1[i] = maxthrust;

        *vecLow->u2[i] = -maxthetadot;
        *vecUpp->u2[i] = maxthetadot;
    }

    //Set boundaries for variables with dimension N-1
    for(int j = 0; j<N_-1; j++){
        *vecLow->u1m[j] = minthrust;
        *vecUpp->u1m[j] = maxthrust;
        *vecLow->u2m[j] = -maxthetadot;
        *vecUpp->u2m[j] = maxthetadot;
    }

    //Constraint vectors for lower and upper boundaries of constraints
    ConstraintVector * gLow = new ConstraintVector(N_);
    ConstraintVector * gUpp = new ConstraintVector(N_);

    //Get the addresses for constraints
    getConstraintVector(gLow,g_l,N_);
    getConstraintVector(gUpp,g_u,N_);

    //Set boundaries for boundary condition constraints and final conditions constraints
    for(int i = 0; i<5 ; i++){
        *gLow->boundCond[i] = 0.;
        *gUpp->boundCond[i] = 0.;
        *gLow->finalCond[i] = 0.;
        *gUpp->finalCond[i] = 0.;
    }

    //Set boundaries for dynamic constraints
    for(int j = 0; j < N_-1; j++){
        *gLow->dynamicX[j] = 0.;
        *gUpp->dynamicX[j] = 0.;

        *gLow->dynamicVX[j] = 0.;
        *gUpp->dynamicVX[j] = 0.;

        *gLow->dynamicZ[j] = 0.;
        *gUpp->dynamicZ[j] = 0.;

        *gLow->dynamicVZ[j] = 0.;
        *gUpp->dynamicVZ[j] = 0.;

        *gLow->dynamicTheta[j] = 0.;
        *gUpp->dynamicTheta[j] = 0.;
    }

//    g_l[m-2] = 0.;
//    g_u[m-2] = 0.;

//    g_l[m-1] = 0.;
//    g_u[m-1] = 0.;

    delete vecLow;
    delete vecUpp;
    delete gLow;
    delete gUpp;

    return true;
}

// returns the initial point for the problem
bool quadNLP::get_starting_point(Index n, bool init_x, Number* x,
                                 bool init_z, Number* z_L, Number* z_U,
                                 Index m, bool init_lambda,
                                 Number* lambda)
{
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);

    //Let tf = tn
    x[9*N_-2] = tn;
    x[0] = this->x0;
    x[2*N_] = this->z0;
    x[N_] = this->vx0;
    x[3*N_] = this->vz0;
    x[4*N_] = this->theta0;
    x[5*N_] = this->u10;
    x[6*N_] = this->u20;

#ifdef skip_me
    /* If checking derivatives, if is useful to  choose different values */
#endif

    return true;
}
// returns the value of the objective function
bool quadNLP::eval_f(Index n, const Number* x,
                     bool new_x, Number& obj_value)
{
    obj_value = x[9*N_-2];
    return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool quadNLP::eval_grad_f(Index n, const Number* x,
                          bool new_x, Number* grad_f)
{
    for(Index i = 0 ; i<n-1; i ++){
        grad_f[i] = 0.;
    }

    grad_f[n-1] = 1.;

    return true;
}

// return the value of the constraints: g(x)
bool quadNLP::eval_g(Index n, const Number* x,
                     bool new_x, Index m, Number* g)
{
    StateVectorConst * sV = new StateVectorConst(N_);
    getStateVectorConst(sV, x, N_);

    ConstraintVector * cV = new ConstraintVector(N_);
    getConstraintVector(cV, g, N_);

    *cV->boundCond[0] = sV->x[0] - this->x0;
    *cV->boundCond[1] = sV->z[0] - this->z0;
    *cV->boundCond[2] = sV->vx[0] - this->vx0;
    *cV->boundCond[3] = sV->vz[0] - this->vz0;
    *cV->boundCond[4] = sV->theta[0] - this->theta0;

    *cV->finalCond[0] = sV->x[N_-1] - xn;
    *cV->finalCond[1] = sV->z[N_-1] - zn;
    *cV->finalCond[2] = sV->vz[N_-1] - vxn;
    *cV->finalCond[3] = sV->vx[N_-1] - vzn;
    *cV->finalCond[4] = sV->theta[N_-1] - thetan;

    Number coeff = sV->tf/(N_-1)/6;

    for(int j = 0 ; j<N_-1 ; j++){
        *cV->dynamicX[j]    =   -sV->x[j]       + computeDynamicX(  sV->x[j+1],     sV->vx[j], sV->vx[j+1],     sV->vxm[j], coeff);
        *cV->dynamicVX[j]   =   -sV->vx[j]      + computeDynamicVX( sV->vx[j+1],    sV->u1[j], sV->theta[j],    sV->u1[j+1], sV->theta[j+1], sV->u1m[j], sV->thetam[j], coeff);
        *cV->dynamicZ[j]    =   -sV->z[j]       + computeDynamicZ(  sV->z[j+1],     sV->vz[j], sV->vz[j+1],     sV->vzm[j], coeff);
        *cV->dynamicVZ[j]   =   -sV->vz[j]      + computeDynamicVZ( sV->vz[j+1],    sV->u1[j], sV->theta[j],    sV->u1[j+1], sV->theta[j+1], sV->u1m[j], sV->thetam[j], coeff);
        *cV->dynamicTheta[j] =  -sV->theta[j]   + computeTheta(     sV->theta[j+1], sV->u2[j], sV->u2[j+1],     sV->u2m[j], coeff);
    }

//    g[m-2] = sV->u1[0] - this->u10;
//    g[m-1] = sV->u2[0] - this->u20;

    delete sV;
    delete cV;
    return true;
}
//nie dziala
// return the structure or values of the jacobian
bool quadNLP::eval_jac_g(Index n, const Number* x, bool new_x,
                         Index m, Index nele_jac, Index* iRow,
                         Index *jCol, Number* values)
{

    Index vx = N_;
    Index z = 2*N_;
    Index vz = 3*N_;
    Index theta = 4*N_;
    Index u1 = 5*N_;
    Index u2 = 6*N_;
    Index u1m = 7*N_;
    Index u2m = 8*N_-1;
    Index tf = 9*N_-2;

    if (values == NULL) {

        // return the structure of the jacobian
        Index ind = 0;
        // X0
        iRow[ind] = 0;
        jCol[ind] = 0;
        ind++;
        // Z0
        iRow[ind] = 1;
        jCol[ind] = z;
        ind++;
        // VX0
        iRow[ind] = 2;
        jCol[ind] = vx;
        ind++;
        // VZ0
        iRow[ind] = 3;
        jCol[ind] = vz;
        ind++;
        // Theta
        iRow[ind] = 4;
        jCol[ind] = theta;
        ind++;
        // XN
        iRow[ind] = 5;
        jCol[ind] = N_-1;
        ind++;
        // ZN
        iRow[ind] = 6;
        jCol[ind] = z+N_-1;
        ind++;
        // VZN
        iRow[ind] = 7;
        jCol[ind] = vz+N_-1;
        ind++;
        // VXN
        iRow[ind] = 8;
        jCol[ind] = vx+N_-1;
        ind++;
        // ThetaN
        iRow[ind] = 9;
        jCol[ind] = theta+N_-1;
        ind++;

        Index lowerRowBound = 10;
        Index upperRowBound = N_-1 + 10;

        Index realInd = 0;

        //DYNAMIC X
        for ( Index j = lowerRowBound; j < upperRowBound; j++){
            realInd = j - lowerRowBound;

            vx = N_ + realInd;
            z = 2*N_+ realInd;
            vz = 3*N_+ realInd;
            theta = 4*N_+ realInd;
            u1 = 5*N_+ realInd;
            u2 = 6*N_+ realInd;
            u1m = 7*N_+ realInd;
            u2m = 8*N_-1 + realInd;

            // X
            iRow[ind] = j;
            jCol[ind] = realInd;
            ind++;

            iRow[ind] = j;
            jCol[ind] = realInd + 1;
            ind++;

            // VX
            iRow[ind] = j;
            jCol[ind] = vx;
            ind++;

            iRow[ind] = j;
            jCol[ind] = vx + 1;
            ind++;

            // Theta
            iRow[ind] = j;
            jCol[ind] = theta;
            ind++;

            iRow[ind] = j;
            jCol[ind] = theta + 1;
            ind++;

            // U1
            iRow[ind] = j;
            jCol[ind] = u1;
            ind++;

            iRow[ind] = j;
            jCol[ind] = u1 + 1;
            ind++;

            // TF
            iRow[ind] = j;
            jCol[ind] = tf;
            ind++;
        }

        lowerRowBound = upperRowBound;
        upperRowBound += N_-1;

        // DYNAMIC VX
        for ( Index j = lowerRowBound; j < upperRowBound; j++){
            realInd = j - lowerRowBound;
            vx = N_ + realInd;
            z = 2*N_+ realInd;
            vz = 3*N_+ realInd;
            theta = 4*N_+ realInd;
            u1 = 5*N_+ realInd;
            u2 = 6*N_+ realInd;
            u1m = 7*N_+ realInd;
            u2m = 8*N_-1+ realInd;

            // VX
            iRow[ind] = j;
            jCol[ind] = vx;
            ind++;

            iRow[ind] = j;
            jCol[ind] = vx + 1;
            ind++;

            // Theta
            iRow[ind] = j;
            jCol[ind] = theta;
            ind++;

            iRow[ind] = j;
            jCol[ind] = theta + 1;
            ind++;

            // U1
            iRow[ind] = j;
            jCol[ind] = u1;
            ind++;

            iRow[ind] = j;
            jCol[ind] = u1 + 1;
            ind++;

            // U2
            iRow[ind] = j;
            jCol[ind] = u2;
            ind++;

            iRow[ind] = j;
            jCol[ind] = u2 + 1;
            ind++;

            // U1M
            iRow[ind] = j;
            jCol[ind] = u1m;
            ind++;

            // TF
            iRow[ind] = j;
            jCol[ind] = tf;
            ind++;
        }

        lowerRowBound = upperRowBound;
        upperRowBound += N_-1;

        // DYNAMIC Z
        for (Index j = lowerRowBound; j < upperRowBound; j++){

            realInd = j - lowerRowBound;
            vx = N_ + realInd;
            z = 2*N_+ realInd;
            vz = 3*N_+ realInd;
            theta = 4*N_+ realInd;
            u1 = 5*N_+ realInd;
            u2 = 6*N_+ realInd;
            u1m = 7*N_+ realInd;
            u2m = 8*N_-1+ realInd;

            // Z
            iRow[ind] = j;
            jCol[ind] = z;
            ind++;

            iRow[ind] = j;
            jCol[ind] = z + 1;
            ind++;

            // VZ
            iRow[ind] = j;
            jCol[ind] = vz;
            ind++;

            iRow[ind] = j;
            jCol[ind] = vz + 1;
            ind++;

            // Theta
            iRow[ind] = j;
            jCol[ind] = theta;
            ind++;

            iRow[ind] = j;
            jCol[ind] = theta + 1;
            ind++;

            // U1
            iRow[ind] = j;
            jCol[ind] = u1;
            ind++;

            iRow[ind] = j;
            jCol[ind] = u1 + 1;
            ind++;

            // TF
            iRow[ind] = j;
            jCol[ind] = tf;
            ind++;
        }
        lowerRowBound = upperRowBound;
        upperRowBound += N_-1;

        //DYNAMIC VZ
        for ( Index j = lowerRowBound; j < upperRowBound; j++){

            realInd = j - lowerRowBound;
            vx = N_ + realInd;
            z = 2*N_+ realInd;
            vz = 3*N_+ realInd;
            theta = 4*N_+ realInd;
            u1 = 5*N_+ realInd;
            u2 = 6*N_+ realInd;
            u1m = 7*N_+ realInd;
            u2m = 8*N_-1+ realInd;

            // VZ
            iRow[ind] = j;
            jCol[ind] = vz;
            ind++;

            iRow[ind] = j;
            jCol[ind] = vz + 1;
            ind++;

            // Theta
            iRow[ind] = j;
            jCol[ind] = theta;
            ind++;

            iRow[ind] = j;
            jCol[ind] = theta + 1;
            ind++;

            // U1
            iRow[ind] = j;
            jCol[ind] = u1;
            ind++;

            iRow[ind] = j;
            jCol[ind] = u1 + 1;
            ind++;

            // U2
            iRow[ind] = j;
            jCol[ind] = u2;
            ind++;

            iRow[ind] = j;
            jCol[ind] = u2 + 1;
            ind++;

            // U1M
            iRow[ind] = j;
            jCol[ind] = u1m;
            ind++;

            // TF
            iRow[ind] = j;
            jCol[ind] = tf;
            ind++;
        }
        lowerRowBound = upperRowBound;
        upperRowBound += N_-1;

        //DYNAMIC THETA
        for ( Index j = lowerRowBound; j < upperRowBound; j++){

            realInd = j - lowerRowBound;
            vx = N_ + realInd;
            z = 2*N_+ realInd;
            vz = 3*N_+ realInd;
            theta = 4*N_+ realInd;
            u1 = 5*N_+ realInd;
            u2 = 6*N_+ realInd;
            u1m = 7*N_+ realInd;
            u2m = 8*N_-1+ realInd;

            // Theta
            iRow[ind] = j;
            jCol[ind] = theta;
            ind++;

            iRow[ind] = j;
            jCol[ind] = theta + 1;
            ind++;

            // U2
            iRow[ind] = j;
            jCol[ind] = u2;
            ind++;

            iRow[ind] = j;
            jCol[ind] = u2 + 1;
            ind++;

            // U2M
            iRow[ind] = j;
            jCol[ind] = u2m;
            ind++;

            // TF
            iRow[ind] = j;
            jCol[ind] = tf;
            ind++;
        }
        lowerRowBound = upperRowBound;
        upperRowBound += N_-1;

//        iRow[ind] = lowerRowBound;
//        jCol[ind] = 5*N_;
//        ind++;

//        iRow[ind] = lowerRowBound+1;
//        jCol[ind] = 6*N_;
//        ind++;

    }
    else {
        // return the values of the jacobian of the constraints

        // assert(inz==nele_jac);
    }

    return true;
}

//return the structure or values of the hessian
bool quadNLP::eval_h(Index n, const Number* x, bool new_x,
                     Number obj_factor, Index m, const Number* lambda,
                     bool new_lambda, Index nele_hess, Index* iRow,
                     Index* jCol, Number* values)
{
    return true;
}

void quadNLP::finalize_solution(SolverReturn status,
                                Index n, const Number* x,
                                const Number* z_L, const Number* z_U,
                                Index m, const Number* g,
                                const Number* lambda,
                                Number obj_value,
                                const IpoptData* ip_data,
                                IpoptCalculatedQuantities* ip_cq)
{
    StateVectorConst * sV = new StateVectorConst(N_);
    getStateVectorConst(sV, x, N_);

    printf("\nWriting solution file solution.txt\n");
    FILE* fp = fopen("solution.txt", "w");

    Number dt = sV->tf/(N_-1);
    Number* timegrid = new Number[N_];
    Number ind = 1;

    for(Index i = 0 ; i <N_; i ++){
        timegrid[i] = dt*(ind-1);
        ind++;
    }

    for (Index i=0; i< N_-1; i++) {
        fprintf(fp, "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", timegrid[i], sV->x[i], sV->vx[i], sV->z[i], sV->vz[i], sV->theta[i], sV->u1[i], sV->u2[i]);
        fprintf(fp, "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", timegrid[i] + dt/2, sV->xm[i], sV->vxm[i], sV->zm[i], sV->vzm[i], sV->thetam[i], sV->u1m[i], sV->u2m[i]);
    }
    fprintf (fp, "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", timegrid[N_-1],  sV->x[N_-1], sV->vx[N_-1], sV->z[N_-1], sV->vz[N_-1], sV->theta[N_-1], sV->u1[N_-1], sV->u2[N_-1]);

    delete sV;
    delete [] timegrid;

    fclose(fp);
}
