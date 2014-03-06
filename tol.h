/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#ifndef _TOL_H_
#define _TOL_H_

#include "datastruct.h"

static int tol_0[] =  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
static int comp_tol_0[] =     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

static int tol_1[] =  {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
			  1, 1, 1, 2, 2, 2, 2, 2, 2, 2,
			  2, 2, 2, 2, 2, 2, 2, 3, 3, 3,
			  3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			  3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
			  4, 4, 4, 4, 4, 4, 4, 4, 4, 5,
			  5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
			  5, 5, 5, 5, 5, 5, 5, 6, 6, 6,
			  6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
			  6, 6, 6, 6, 6, 7, 7, 7, 7, 7};

static int comp_tol_1[] =     {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
				  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
				  1, 1, 1, 1, 1, 1, 1, 2, 2, 2,
				  2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
				  2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
				  2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
				  3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
				  3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
				  3, 3, 3, 3, 3, 3, 3, 3, 4, 4,
				  4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

static int tol_2[] =  {0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
			  2, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			  4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
			  4, 4, 5, 5, 5, 5, 5, 5, 5, 5,
			  5, 5, 5, 5, 6, 6, 6, 6, 6, 6,
			  6, 6, 6, 6, 6, 6, 6, 6, 7, 7,
			  7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
			  7, 7, 8, 8, 8, 8, 8, 8, 8, 8,
			  8, 8, 8, 8, 8, 8, 9, 9, 9, 9,
			  9, 9, 9, 9, 9, 9, 9, 9, 9, 9};
static int comp_tol_2[] =     {0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2,
				  2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
				  2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
				  3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
				  3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
				  4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				  4, 4, 4, 4, 4, 5, 5, 5, 5, 5,
				  5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
				  5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
				  5, 6, 6, 6, 6, 6, 6, 6, 6, 6};

int* selectTol(filter_pars_t* pars);
int* selectCompTol(filter_pars_t* pars);
int* loadTol(filter_pars_t* pars);
int* loadCompTol(filter_pars_t* pars);
void printTol(int* tol, int* comp_tol);

#endif