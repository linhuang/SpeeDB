/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#ifndef _DBS_H_
#define _DBS_H_

#include "datastruct.h"

dbs_t* newDBS(datasize_t* data_size);
void freeDBS(dbs_t* dbs);
void storeDBS(char* output_filename, dbs_t* dbs);
void printDBS(dbs_t* dbs);
dbs_t* loadDBS(char* dbs_filename);
void setBlockSNP(dbs_t* dbs, filter_pars_t* pars);
void setStepSNP(dbs_t* dbs, filter_pars_t* pars);

#endif
