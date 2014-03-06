/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#ifndef _QUERY_H_
#define _QUERY_H_

#include "datastruct.h"

queries_t* newQuery(dbs_t* dbs, int query_count);
void printQuery(dbs_t* dbs, queries_t* queries);
void freeQuery(queries_t* queries);
void genQUERY(char* input_filename, char* output_filename);
void genQ(char* input_filename, char* output_filename, query_pars_t* pars);
void storeQuery(char* output_filename, dbs_t* dbs, queries_t* queries);
queries_t* loadQuery(char* query_filename, dbs_t* dbs, filter_pars_t* pars);

#endif
