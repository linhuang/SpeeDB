/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#ifndef _PRECALC_H_
#define _PRECALC_H_

#include "datastruct.h"

precalc_t* newPrecalc(dbs_t* dbs);
void freePrecalc(dbs_t* dbs, precalc_t* precalc);
void printPrecalc(dbs_t* dbs, precalc_t* precalc);
precalc_t* loadPrecalc(dbs_t* dbs, char* precalc_filename);
precalc_t* precalcDBS(dbs_t* dbs, char* precalc_filename, filter_pars_t* pars);

#endif
