/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#ifndef _STEP_H_
#define _STEP_H_

#include "datastruct.h"

steps_t* newSteps(dbs_t* dbs, int64_t steps_count);
void freeSteps(steps_t* steps);
int64_t countSteps(dbs_t* dbs, filter_pars_t* pars);
void setSteps(dbs_t* dbs, steps_t* steps, filter_pars_t* pars);
void printSteps(steps_t* steps, dbs_t* dbs);

#endif
