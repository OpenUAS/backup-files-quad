#ifndef QUADFUNCTIONS_H
#define QUADFUNCTIONS_H
#include <coin/IpTNLP.hpp>

using namespace Ipopt;

struct StateVectorConst{
    Number tf;
    Number* x;
    Number* vx;
    Number* z;
    Number* vz;
    Number* theta;
    Number* u1;
    Number* u2;
    Number* u1m;
    Number* u2m;
    Number* xm;
    Number* vxm;
    Number* zm;
    Number* vzm;
    Number* thetam;
    Index N;

    StateVectorConst(Index N){
        this->N = N;
        x = new Number[N];
        vx = new Number[N];
        z = new Number[N];
        vz = new Number[N];
        theta = new Number[N];
        u1 = new Number[N];
        u2 = new Number[N];

        u1m = new Number[N-1];
        u2m = new Number[N-1];
        xm = new Number[N-1];
        vxm = new Number[N-1];
        zm = new Number[N-1];
        vzm = new Number[N-1];
        thetam = new Number[N-1];
    }

    ~StateVectorConst(){
        delete [] x;
        delete [] vx;
        delete [] z;
        delete [] vz;
        delete [] theta;
        delete [] u1;
        delete [] u2;

        delete [] u1m;
        delete [] u2m;
        delete [] xm;
        delete [] vxm;
        delete [] zm;
        delete [] vzm;
        delete [] thetam;
    }
};

struct StateVector{
    Number* tf;
    Number** x;
    Number** vx;
    Number** z;
    Number** vz;
    Number** theta;
    Number** u1;
    Number** u2;
    Number** u1m;
    Number** u2m;
    Number* xm;
    Number* vxm;
    Number* zm;
    Number* vzm;
    Number* thetam;
    Index N;

    StateVector(Index N){
        this->N = N;
        x = new Number*[N];
        vx = new Number*[N];
        z = new Number*[N];
        vz = new Number*[N];
        theta = new Number*[N];
        u1 = new Number*[N];
        u2 = new Number*[N];

        u1m = new Number*[N-1];
        u2m = new Number*[N-1];

        xm = new Number[N-1];
        vxm = new Number[N-1];
        zm = new Number[N-1];
        vzm = new Number[N-1];
        thetam = new Number[N-1];
    }

    ~StateVector(){
        delete [] x;
        delete [] vx;
        delete [] z;
        delete [] vz;
        delete [] theta;
        delete [] u1;
        delete [] u2;

        delete [] u1m;
        delete [] u2m;

        delete [] xm;
        delete [] vxm;
        delete [] zm;
        delete [] vzm;
        delete [] thetam;
    }
};

struct ConstraintVector{
    Number** boundCond;
    Number** finalCond;

    Number** dynamicX;
    Number** dynamicVX;
    Number** dynamicZ;
    Number** dynamicVZ;
    Number** dynamicTheta;
    Index N;

    ConstraintVector(Index N){
        this->N = N;

        boundCond = new Number*[5];
        finalCond = new Number*[5];

        dynamicX = new Number*[N-1];
        dynamicVX = new Number*[N-1];
        dynamicZ = new Number*[N-1];
        dynamicVZ = new Number*[N-1];
        dynamicTheta = new Number*[N-1];

    }

    ~ConstraintVector(){
        delete [] boundCond;
        delete [] finalCond;

        delete [] dynamicX;
        delete [] dynamicVX;
        delete [] dynamicZ;
        delete [] dynamicVZ;
        delete [] dynamicTheta;
    }
};

Number computeDynamicX(Number,Number,Number,Number,Number);
Number computeDynamicVX(Number,Number,Number,Number,Number,Number,Number,Number);
Number computeDynamicZ(Number,Number,Number,Number,Number);
Number computeDynamicVZ(Number,Number,Number,Number,Number,Number,Number,Number);
Number computeTheta(Number,Number,Number,Number,Number);

Number computeXM(Number,Number,Number,Number,Number);
Number computeVXM(Number,Number,Number,Number,Number,Number,Number);
Number computeZM(Number,Number,Number,Number,Number);
Number computeVZM(Number,Number,Number,Number,Number,Number);
Number computeThetaM(Number,Number,Number,Number,Number);

void getStateVectorConst(StateVectorConst *,const Number *, Index);
void getStateVector(StateVector *, Number *, Index);
void getConstraintVector(ConstraintVector *,Number *, Index);

#endif // QUADFUNCTIONS_H
