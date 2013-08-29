/*
 * optimal_trajectory.h
 *
 *  Created on: Aug 21, 2013
 *      Author: bheliom
 */

#ifndef OPTIMAL_TRAJECTORY_H_
#define OPTIMAL_TRAJECTORY_H_

#include "IpStdCInterface.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

/* Function Declarations */
Bool eval_f(Index n, Number* x, Bool new_x, Number* obj_value,
		UserDataPtr user_data);

Bool eval_grad_f(Index n, Number* x, Bool new_x, Number* grad_f,
		UserDataPtr user_data);

Bool eval_g(Index n, Number* x, Bool new_x, Index m, Number* g,
		UserDataPtr user_data);

Bool eval_jac_g(Index n, Number *x, Bool new_x, Index m, Index nele_jac,
		Index *iRow, Index *jCol, Number *values, UserDataPtr user_data);

Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor, Index m,
		Number *lambda, Bool new_lambda, Index nele_hess, Index *iRow,
		Index *jCol, Number *values, UserDataPtr user_data);
int ipoptMain(void);

extern float ipopt_cmd[2];

#endif /* OPTIMAL_TRAJECTORY_H_ */
