/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#ifndef _FILT_H_
#define _FILT_H_

#include "datastruct.h"

precalc_t* precalculation(char* dbs_filename, dbs_t* dbs, filter_pars_t* pars);
int filterDB(char* dbs_filename, char* query_filename, char* output_filename, filter_pars_t* pars);

#endif
