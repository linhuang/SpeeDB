/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#ifndef _COMMON_H_
#define _COMMON_H_

#include <stdio.h>
#include <stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#define ZERO 1e-8

int isMissing(char c);
double getMaxDouble(double a, double b);
double getMinDouble(double a, double b);
int64_t getMaxInt64(int64_t a, int64_t b);
FILE* fileOpenR(char* filename);
FILE* fileOpenRB(char* filename);
FILE* fileOpenW(char* filename);
FILE* fileOpenWB(char* filename);
void fileClose(FILE* file);
int fileExists(const char *filename);

#endif
