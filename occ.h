/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#ifndef _OCC_H_
#define _OCC_H_

#include "datastruct.h"

occ_t* newOcc(dbs_t* dbs);
void clearOcc(occ_t* occ, dbs_t* dbs);
void printOcc(occ_t* occ, dbs_t* dbs);
void freeOcc(occ_t* occ);

#endif
