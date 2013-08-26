/*
 * quadNLPmain.c
 *
 *  Created on: Aug 20, 2013
 *      Author: bheliom
 */

#include "optimal_trajectory.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <state.h>
#include <math/pprz_algebra_int.h>

float ipopt_cmd[2];

/*Number of nodes*/
Index N_ = 7;

/*Starting indeces in x array for variables*/
Index xP;
Index vxP;
Index zP;
Index vzP;
Index thetaP;
Index u1P;
Index u2P;
Index u1mP;
Index u2mP;

Number x0 = 0;
Number z0 = 0;
Number vx0 = 0.0;
Number vz0 = 0.0;
Number theta0 = 0.0;

Number xn = 25.0;
Number zn = 10.0;
Number vxn = 0.0;
Number vzn = 0.0;
Number thetan = 0.0;

Number maxthrust = 20.0;
Number minthrust = 1.0;
Number maxthetadot = 5.0;
Number m0 = 1.0;
Number gearth = 9.81;

Number tn = 1.0;

struct MyUserData {
	Number g_offset[2]; /* This is an offset for the constraints.  */
};

int ipoptMain(void) {

	Index n = -1; /* number of variables */
	Index m = -1; /* number of constraints */
	Number* x_L = NULL; /* lower bounds on x */
	Number* x_U = NULL; /* upper bounds on x */
	Number* g_L = NULL; /* lower bounds on g */
	Number* g_U = NULL; /* upper bounds on g */
	IpoptProblem nlp = NULL; /* IpoptProblem */
	Index upperB;

	enum ApplicationReturnStatus status; /* Solve return code */

	Number* x = NULL; /* starting point and solution vector */
	Number* mult_g = NULL; /* constraint multipliers
	 at the solution */

	Number* mult_x_L = NULL; /* lower bound multipliers
	 at the solution */
	Number* mult_x_U = NULL; /* upper bound multipliers
	 at the solution */

	Number obj; /* objective value */

	Index i;
	Index j;

	struct MyUserData user_data;
	/* Number of nonzeros in the Jacobian of the constraints */
	Index nele_jac = 44 * (N_ - 1) + 10;

	Index nele_hess = 0;
	/* indexing style for matrices */
	Index index_style = 0; /* C-style; start counting of rows and column
	 indices at 0 */

	/* set the number of variables and allocate space for the bounds */
	n = 9 * N_ - 1;
	x_L = (Number*) malloc(sizeof(Number) * n);
	x_U = (Number*) malloc(sizeof(Number) * n);

	x_L[n - 1] = 0.0;
	x_U[n - 1] = 2e19;

	xP = 0;
	vxP = N_;
	zP = 2 * N_;
	vzP = 3 * N_;
	thetaP = 4 * N_;
	u1P = 5 * N_;
	u2P = 6 * N_;
	u1mP = 7 * N_;
	u2mP = 8 * N_ - 1;

	for (i = 0; i < N_; i++) {
		x_L[xP + i] = -2e19;
		x_U[xP + i] = 2e19;

		x_L[vxP + i] = -2e19;
		x_U[vxP + i] = 2e19;

		x_L[zP + i] = -2e19;
		x_U[zP + i] = 2e19;

		x_L[vzP + i] = -2e19;
		x_U[vzP + i] = 2e19;

		x_L[thetaP + i] = -2e19;
		x_U[thetaP + i] = 2e19;

		x_L[u1P + i] = minthrust;
		x_U[u1P + i] = maxthrust;

		x_L[u2P + i] = -maxthetadot;
		x_U[u2P + i] = maxthetadot;
	}

	for (j = 0; j < N_ - 1; j++) {
		x_L[u1mP + j] = minthrust;
		x_U[u1mP + j] = maxthrust;
		x_L[u2mP + j] = -maxthetadot;
		x_U[u2mP + j] = maxthetadot;
	}

	/* set the number of constraints and allocate space for the bounds */
	m = 5 * (N_ - 1) + 10;
	g_L = (Number*) malloc(sizeof(Number) * m);
	g_U = (Number*) malloc(sizeof(Number) * m);

	for (i = 0; i < 5; i++) {
		g_L[i] = 0.0;
		g_U[i] = 0.0;
		g_L[i + 5] = 0.0;
		g_U[i + 5] = 0.0;
	}

	upperB = N_ - 1;

	for (j = 0; j < upperB; j++) {
		g_L[10 + j] = 0.0;
		g_U[10 + j] = 0.0;

		g_L[10 + upperB + j] = 0.0;
		g_U[10 + upperB + j] = 0.0;

		g_L[10 + 2 * upperB + j] = 0.0;
		g_U[10 + 2 * upperB + j] = 0.0;

		g_L[10 + 3 * upperB + j] = 0.0;
		g_U[10 + 3 * upperB + j] = 0.0;

		g_L[10 + 4 * upperB + j] = 0.0;
		g_U[10 + 4 * upperB + j] = 0.0;
	}

	/* create the IpoptProblem */

	nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
			index_style, &eval_f, &eval_g, &eval_grad_f, &eval_jac_g, &eval_h);

	/* We can free the memory now - the values for the bounds have been
	 copied internally in CreateIpoptProblem */
	free(x_L);
	free(g_L);
	free(x_U);
	free(g_U);

	/* Set some options.  Note the following ones are only examples,
	 they might not be suitable for your problem. */

	AddIpoptNumOption(nlp, "tol", 1e-3);
	AddIpoptStrOption(nlp, "jacobian_approximation",
			"finite-difference-values");
	AddIpoptIntOption(nlp, "max_iter", 1000);
	AddIpoptStrOption(nlp, "hessian_approximation", "limited-memory");
	AddIpoptIntOption(nlp, "print_level", 3);

	/* allocate space for the initial point and set the values */
	x = (Number*) malloc(sizeof(Number) * n);
	x[n - 1] = tn;

	for (i = 0; i < N_; i++) {
		x[xP + i] = x0;
		x[zP + i] = z0;
		x[vxP + i] = vx0;
		x[vzP + i] = vz0;
		x[thetaP + i] = theta0;
	}
	/* allocate space to store the bound multipliers at the solution */
	mult_g = (Number*) malloc(sizeof(Number) * m);
	mult_x_L = (Number*) malloc(sizeof(Number) * n);
	mult_x_U = (Number*) malloc(sizeof(Number) * n);

	/* Initialize the user data */
	user_data.g_offset[0] = 0.;
	user_data.g_offset[1] = 0.;

	/* Set the callback method for intermediate user-control.  This is
	 * not required, just gives you some intermediate control in case
	 * you need it. */
	/* SetIntermediateCallback(nlp, intermediate_cb); */
	/* solve the problem */
	status = IpoptSolve(nlp, x, NULL, &obj, mult_g, mult_x_L, mult_x_U,
			&user_data);

	ipopt_cmd[0] = x[u1P+1];
	ipopt_cmd[1] = x[thetaP+1] - x[thetaP];

	return (int) status;
}

/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
/*------------------------------------------FUNCTIONS--------------------------------------------*/
Number computeDynamicX(Number x1, Number vx, Number vx1, Number vxm,
		Number coeff) {
	return x1 - coeff * (vx + vx1 + 4 * vxm);
}

Number computeDynamicVX(Number vx1, Number u1, Number theta, Number u1_1,
		Number theta_1, Number u1m, Number thetam, Number coeff) {
	return vx1
			- coeff
					* (u1 * sin(theta) / m0 + u1_1 * sin(theta_1) / m0
							+ 4 * (u1m * sin(thetam) / m0));
}

Number computeDynamicZ(Number z1, Number vz, Number vz1, Number vzm,
		Number coeff) {
	return z1 - coeff * (vz + vz1 + 4 * vzm);
}

Number computeDynamicVZ(Number vz1, Number u1, Number theta, Number u1_1,
		Number theta1, Number u1m, Number thetam, Number coeff) {
	return vz1
			- coeff
					* ((u1 * cos(theta) / m0 - gearth)
							+ (u1_1 * cos(theta1) / m0 - gearth)
							+ 4 * (u1m * cos(thetam) / m0 - gearth));
}

Number computeTheta(Number theta1, Number u2, Number u2_1, Number u2m,
		Number coeff) {
	return theta1 - coeff * (u2 + u2_1 + 4 * u2m);
}

Number computeVXM(Number vx, Number vx1, Number u1, Number theta, Number u1_1,
		Number theta1, Number coeff) {
	return (vx + vx1) / 2
			+ coeff * (u1 * sin(theta) / m0 - u1_1 * sin(theta1) / m0);
}

Number computeVZM(Number vz, Number vz1, Number u1, Number theta, Number u1_1,
		Number theta1, Number coeff) {
	return (vz + vz1) / 2
			+ coeff
					* ((u1 * cos(theta) / m0 - gearth)
							- (u1_1 * cos(theta1) / m0 - gearth));
}

Number computeThetaM(Number theta, Number theta1, Number u2, Number u2_1,
		Number coeff) {
	return (theta + theta1) / 2 + coeff * (u2 - u2_1);
}

Bool eval_f(Index n, Number* x, Bool new_x, Number* obj_value,
		UserDataPtr user_data) {

	*obj_value = x[n - 1];

	return TRUE;
}

Bool eval_grad_f(Index n, Number* x, Bool new_x, Number* grad_f,
		UserDataPtr user_data) {
	Index i;
	for (i = 0; i < n - 1; i++) {
		grad_f[i] = 0.0;
	}

	grad_f[n - 1] = 1.0;

	return TRUE;
}

Bool eval_g(Index n, Number* x, Bool new_x, Index m, Number* g,
		UserDataPtr user_data) {

	/*struct MyUserData* my_data = user_data;*/
	Number coeff1;
	Number coeff2;
	Index upperB;
	Index j;

	struct FloatEulers* currAng = stateGetNedToBodyEulers_f();
	struct EnuCoor_f* currCord = stateGetPositionEnu_f();
	struct EnuCoor_f* currSpeed = stateGetSpeedEnu_f();

	g[0] = x[xP] - currCord->x;
	g[1] = x[zP] - currCord->z;
	g[2] = x[vxP] - currSpeed->x;
	g[3] = x[vzP] - currSpeed->z;
	g[4] = x[thetaP] - currAng->theta;

	g[5] = x[xP + N_ - 1] - xn;
	g[6] = x[zP + N_ - 1] - zn;
	g[7] = x[vzP + N_ - 1] - vzn;
	g[8] = x[vxP + N_ - 1] - vxn;
	g[9] = x[thetaP + N_ - 1] - thetan;

	coeff1 = x[n - 1] / (N_ - 1) / 6;
	coeff2 = x[n - 1] / (N_ - 1) / 8;
	upperB = N_ - 1;

	for (j = 0; j < upperB; j++) {
		g[10 + j] = -x[xP + j]
				+ computeDynamicX(x[xP + j + 1], x[vxP + j], x[vxP + j + 1],
						computeVXM(x[vxP + j], x[vxP + j + 1], x[u1P + j],
								x[thetaP + j], x[u1P + j + 1],
								x[thetaP + j + 1], coeff2), coeff1);

		g[10 + upperB + j] = -x[vxP + j]
				+ computeDynamicVX(x[vxP + j + 1], x[u1P + j], x[thetaP + j],
						x[u1P + j + 1], x[thetaP + j + 1], x[u1mP + j],
						computeThetaM(x[thetaP + j], x[thetaP + j + 1],
								x[u2P + j], x[u2P + j + 1], coeff2), coeff1);
		g[10 + 2 * upperB + j] = -x[zP + j]
				+ computeDynamicZ(x[zP + j + 1], x[vzP + j], x[vzP + j + 1],
						computeVZM(x[vzP + j], x[vzP + j + 1], x[u1P + j],
								x[thetaP + j], x[u1P + j + 1],
								x[thetaP + j + 1], coeff2), coeff1);
		g[10 + 3 * upperB + j] = -x[vzP + j]
				+ computeDynamicVZ(x[vzP + j + 1], x[u1P + j], x[thetaP + j],
						x[u1P + j + 1], x[thetaP + j + 1], x[u1mP + j],
						computeThetaM(x[thetaP + j], x[thetaP + j + 1],
								x[u2P + j], x[u2P + j + 1], coeff2), coeff1);
		g[10 + 4 * upperB + j] = -x[thetaP + j]
				+ computeTheta(x[thetaP + j + 1], x[u2P + j], x[u2P + j + 1],
						x[u2mP + j], coeff1);
	}

	return TRUE;
}

Bool eval_jac_g(Index n, Number *x, Bool new_x, Index m, Index nele_jac,
		Index *iRow, Index *jCol, Number *values, UserDataPtr user_data) {

	Index vx = N_;
	Index z = 2 * N_;
	Index vz = 3 * N_;
	Index theta = 4 * N_;
	Index u1 = 5 * N_;
	Index u2 = 6 * N_;
	Index u1m = 7 * N_;
	Index u2m = 8 * N_ - 1;
	Index tf = 9 * N_ - 2;
	Index lowerRowBound;
	Index upperRowBound;
	Index realInd;
	Index j;
	Index ind;

	if (values == NULL ) {
		ind = 0;

		iRow[ind] = 0;
		jCol[ind] = 0;
		ind++;

		iRow[ind] = 1;
		jCol[ind] = z;
		ind++;

		iRow[ind] = 2;
		jCol[ind] = vx;
		ind++;

		iRow[ind] = 3;
		jCol[ind] = vz;
		ind++;

		iRow[ind] = 4;
		jCol[ind] = theta;
		ind++;

		iRow[ind] = 5;
		jCol[ind] = N_ - 1;
		ind++;

		iRow[ind] = 6;
		jCol[ind] = z + N_ - 1;
		ind++;

		iRow[ind] = 7;
		jCol[ind] = vz + N_ - 1;
		ind++;

		iRow[ind] = 8;
		jCol[ind] = vx + N_ - 1;
		ind++;

		iRow[ind] = 9;
		jCol[ind] = theta + N_ - 1;
		ind++;

		lowerRowBound = 10;
		upperRowBound = N_ - 1 + 10;

		realInd = 0;

		for (j = lowerRowBound; j < upperRowBound; j++) {
			realInd = j - lowerRowBound;

			vx = N_ + realInd;
			z = 2 * N_ + realInd;
			vz = 3 * N_ + realInd;
			theta = 4 * N_ + realInd;
			u1 = 5 * N_ + realInd;
			u2 = 6 * N_ + realInd;
			u1m = 7 * N_ + realInd;
			u2m = 8 * N_ - 1 + realInd;

			iRow[ind] = j;
			jCol[ind] = realInd;
			ind++;

			iRow[ind] = j;
			jCol[ind] = realInd + 1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = vx;
			ind++;

			iRow[ind] = j;
			jCol[ind] = vx + 1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = theta;
			ind++;

			iRow[ind] = j;
			jCol[ind] = theta + 1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = u1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = u1 + 1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = tf;
			ind++;
		}

		lowerRowBound = upperRowBound;
		upperRowBound += N_ - 1;

		for (j = lowerRowBound; j < upperRowBound; j++) {
			realInd = j - lowerRowBound;
			vx = N_ + realInd;
			z = 2 * N_ + realInd;
			vz = 3 * N_ + realInd;
			theta = 4 * N_ + realInd;
			u1 = 5 * N_ + realInd;
			u2 = 6 * N_ + realInd;
			u1m = 7 * N_ + realInd;
			u2m = 8 * N_ - 1 + realInd;

			iRow[ind] = j;
			jCol[ind] = vx;
			ind++;

			iRow[ind] = j;
			jCol[ind] = vx + 1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = theta;
			ind++;

			iRow[ind] = j;
			jCol[ind] = theta + 1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = u1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = u1 + 1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = u2;
			ind++;

			iRow[ind] = j;
			jCol[ind] = u2 + 1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = u1m;
			ind++;

			iRow[ind] = j;
			jCol[ind] = tf;
			ind++;
		}

		lowerRowBound = upperRowBound;
		upperRowBound += N_ - 1;

		for (j = lowerRowBound; j < upperRowBound; j++) {

			realInd = j - lowerRowBound;
			vx = N_ + realInd;
			z = 2 * N_ + realInd;
			vz = 3 * N_ + realInd;
			theta = 4 * N_ + realInd;
			u1 = 5 * N_ + realInd;
			u2 = 6 * N_ + realInd;
			u1m = 7 * N_ + realInd;
			u2m = 8 * N_ - 1 + realInd;

			iRow[ind] = j;
			jCol[ind] = z;
			ind++;

			iRow[ind] = j;
			jCol[ind] = z + 1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = vz;
			ind++;

			iRow[ind] = j;
			jCol[ind] = vz + 1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = theta;
			ind++;

			iRow[ind] = j;
			jCol[ind] = theta + 1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = u1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = u1 + 1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = tf;
			ind++;
		}
		lowerRowBound = upperRowBound;
		upperRowBound += N_ - 1;

		for (j = lowerRowBound; j < upperRowBound; j++) {

			realInd = j - lowerRowBound;
			vx = N_ + realInd;
			z = 2 * N_ + realInd;
			vz = 3 * N_ + realInd;
			theta = 4 * N_ + realInd;
			u1 = 5 * N_ + realInd;
			u2 = 6 * N_ + realInd;
			u1m = 7 * N_ + realInd;
			u2m = 8 * N_ - 1 + realInd;

			iRow[ind] = j;
			jCol[ind] = vz;
			ind++;

			iRow[ind] = j;
			jCol[ind] = vz + 1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = theta;
			ind++;

			iRow[ind] = j;
			jCol[ind] = theta + 1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = u1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = u1 + 1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = u2;
			ind++;

			iRow[ind] = j;
			jCol[ind] = u2 + 1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = u1m;
			ind++;

			iRow[ind] = j;
			jCol[ind] = tf;
			ind++;
		}
		lowerRowBound = upperRowBound;
		upperRowBound += N_ - 1;

		for (j = lowerRowBound; j < upperRowBound; j++) {

			realInd = j - lowerRowBound;
			vx = N_ + realInd;
			z = 2 * N_ + realInd;
			vz = 3 * N_ + realInd;
			theta = 4 * N_ + realInd;
			u1 = 5 * N_ + realInd;
			u2 = 6 * N_ + realInd;
			u1m = 7 * N_ + realInd;
			u2m = 8 * N_ - 1 + realInd;

			iRow[ind] = j;
			jCol[ind] = theta;
			ind++;

			iRow[ind] = j;
			jCol[ind] = theta + 1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = u2;
			ind++;

			iRow[ind] = j;
			jCol[ind] = u2 + 1;
			ind++;

			iRow[ind] = j;
			jCol[ind] = u2m;
			ind++;

			iRow[ind] = j;
			jCol[ind] = tf;
			ind++;
		}
	}
	return TRUE;
}

Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor, Index m,
		Number *lambda, Bool new_lambda, Index nele_hess, Index *iRow,
		Index *jCol, Number *values, UserDataPtr user_data) {
	return TRUE;
}

