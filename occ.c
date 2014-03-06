/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "occ.h"

occ_t* newOcc(dbs_t* dbs)
{
        occ_t* occ = (occ_t*) malloc(1 * sizeof(occ_t));
        if(occ == NULL) printf("malloc error: occ\n");
        if((occ->occ_fw = (int*) malloc(dbs->tuple_count * sizeof(int))) == NULL) printf("malloc error: occ\n");
        if((occ->occ_cp = (int*) malloc(dbs->tuple_count * sizeof(int))) == NULL) printf("malloc error: occ\n");
        return occ;
}

void clearOcc(occ_t* occ, dbs_t* dbs)
{
        memset(occ->occ_fw, 0, dbs->tuple_count * sizeof(int));
        memset(occ->occ_cp, 0, dbs->tuple_count * sizeof(int));
}

void printOcc(occ_t* occ, dbs_t* dbs)
{
        int64_t i;
        for(i = 0; i < dbs->tuple_count; i++)
        {
                printf("%" PRIu64 "(%d), ", i, occ->occ_fw[i]);
        }
        printf("\n");
        for(i = 0; i < dbs->tuple_count; i++)
        {
                printf("%" PRIu64 "(%d), ", i, occ->occ_cp[i]);
        }
        printf("\n");
}

void freeOcc(occ_t* occ)
{
        free(occ->occ_fw);
        free(occ->occ_cp);
        free(occ);
}
