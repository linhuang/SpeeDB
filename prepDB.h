/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#ifndef _PREP_DB_H_
#define _PREP_DB_H_

#include "datastruct.h"

#define TPED_HEADER 4
#define MAX_GEN_POS_LENGTH 20

datasize_t* countSNPs(char* input_filename);
dbs_t* identifyAlleles(char* input_filename, datasize_t* data_size);
int genDB(char* input_filename, char* output_filename);

#endif
