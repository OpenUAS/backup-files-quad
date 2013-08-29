#include <cstdio>

#include <coin/IpIpoptApplication.hpp>
#include "quadNLP.h"
#include <time.h>

using namespace Ipopt;

int main(int argv, char* argc[]){

    // Number of variables
    clock_t sTime = clock();

    Index N = atoi(argc[1]);
    Number x0 = atof(argc[2]);
    Number z0 = atof(argc[3]);
    Number vx0 = atof(argc[4]);
    Number vz0 = atof(argc[5]);
    Number theta0 = atof(argc[6]);

    Number xn = atof(argc[7]);
    Number zn = atof(argc[8]);
    Number thetan = atof(argc[9]);
    Number u10 = atof(argc[10]);
    Number u20 = atof(argc[11]);

    SmartPtr<TNLP> mynlp = new quadNLP(N, x0,z0,vx0,vz0,theta0,xn,zn,thetan,u10,u20);
    SmartPtr<IpoptApplication> app = new IpoptApplication();

    app->Options()->SetNumericValue("tol", 1e-2);
    app->Options()->SetStringValue("jacobian_approximation","finite-difference-values");

    //  app->Options()->SetStringValue("derivative_test","first-order");

    app->Options()->SetIntegerValue("max_iter", 100);
    app->Options()->SetStringValue("hessian_approximation","limited-memory");

    app->Initialize();

    // Ask Ipopt to solve the problem
    ApplicationReturnStatus status = app->OptimizeTNLP(mynlp);

    if (status == Solve_Succeeded) {
        printf("\n\n*** The problem solved!\n");
    }
    else {
        printf("\n\n*** The problem FAILED!\n");
    }

    printf("\nTotal CPU time used: %f", (double) (clock()-sTime)/CLOCKS_PER_SEC," seconds\n\n" );
    return (int) status;
}
