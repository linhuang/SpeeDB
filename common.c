/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#include <stdlib.h>
#include "common.h"

int isMissing(char c)
{
        if(c == 'A' || c == 'C' || c == 'G' || c == 'T') return 0;
        return 1;
}

double getMaxDouble(double a, double b)
{
	if(a > b) return a; else return b;
}

double getMinDouble(double a, double b)
{
	if(a < b) return a; else return b;
}

int64_t getMaxInt64(int64_t a, int64_t b)
{
	if(a > b) return a; else return b;
}

FILE* fileOpenR(char* filename)
{
        FILE* file = (FILE*) fopen(filename, "r");
        if(file == NULL)
        {
                printf("Cannot open the file %s\n", filename);
                exit(1);
        }
        return file;
}

FILE* fileOpenRB(char* filename)
{
        FILE* file = (FILE*) fopen(filename, "rb");
        if(file == NULL)
        {
                printf("Cannot open the file %s\n", filename);
                exit(1);
        }
        return file;
}

FILE* fileOpenW(char* filename)
{
        FILE* file = (FILE*) fopen(filename, "w");
        if(file == NULL)
        {
                printf("Cannot open the file %s\n", filename);
                exit(1);
        }
        return file;
}

FILE* fileOpenWB(char* filename)
{
        FILE* file = (FILE*) fopen(filename, "wb");
        if(file == NULL)
        {
                printf("Cannot open the file %s\n", filename);
                exit(1);
        }
        return file;
}

void fileClose(FILE* file)
{
        fclose(file);
}

int fileExists(const char *filename)
{
        FILE *file = fopen(filename, "r");
        if(file)
        {
                fclose(file);
                return 1;
        }
        return 0;
}
