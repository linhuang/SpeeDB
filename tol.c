/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "common.h"
#include "tol.h"

int* selectTol(filter_pars_t* pars)
{
        if(pars->p_th_fw == 0) return tol_0;
        else if(pars->p_th_fw == 1) return tol_1;
        else if(pars->p_th_fw == 2) return tol_2;
        else return NULL;
}

int* selectCompTol(filter_pars_t* pars)
{
        if(pars->p_th_cp == 0) return comp_tol_0;
        else if(pars->p_th_cp == 1) return comp_tol_1;
        else if(pars->p_th_cp == 2) return comp_tol_2;
        else return NULL;
}

int* loadTol(filter_pars_t* pars)
{
        int* tol = (int*) malloc(101 * sizeof(int));

        char* input_filename = (char*)malloc(strlen(pars->p_filename) + 3);
        sprintf(input_filename, "%s.1", pars->p_filename);
        FILE* input = fileOpenR(input_filename);

        printf("Load filter 1's p threshold from file %s\n", input_filename);
        tol[0] = 0;
        int i, ret;
        for(i = 1; i <= 100; i++)
        {
                ret = fscanf (input, "%d", &tol[i]);
        }

	free(input_filename);
        fileClose(input);

        return tol;
}

int* loadCompTol(filter_pars_t* pars)
{
        int* comp_tol = (int*) malloc(101 * sizeof(int));

        char* input_filename = (char*)malloc(strlen(pars->q_filename) + 3);
        sprintf(input_filename, "%s.2", pars->q_filename);
        FILE* input = fileOpenR(input_filename);

        printf("Load filter 2's p threshold from file %s\n", input_filename);
        comp_tol[0] = 0;
        int i, ret;
        for(i = 1; i <= 100; i++)
        {
                ret = fscanf (input, "%d", &comp_tol[i]);
        }

	free(input_filename);
        fileClose(input);

        return comp_tol;
}

void printTol(int* tol, int* comp_tol)
{
        int i;
        printf("tol:\n");
        for(i = 1; i <= 100; i++)
        {
                printf("%d, ", tol[i]);
                if(i % 10 == 0) printf("\n");
        }
        printf("comp_tol:\n");
        for(i = 1; i <= 100; i++)
        {
                printf("%d, ", comp_tol[i]);
                if(i % 10 == 0) printf("\n");
        }
}
