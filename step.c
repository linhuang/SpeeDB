/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "step.h"

steps_t* newSteps(dbs_t* dbs, int64_t steps_count)
{
        steps_t* steps = (steps_t*) malloc(1 * sizeof(steps_t));
        steps->steps_count = steps_count;
        if((steps->start_genetic_position = (double*) malloc(steps->steps_count * sizeof(double))) == NULL) printf("malloc error: steps\n");
        if((steps->end_genetic_position = (double*) malloc(steps->steps_count * sizeof(double))) == NULL) printf("malloc error: steps\n");
	if((steps->start_SNP = (int64_t*) malloc(steps->steps_count * sizeof(int64_t))) == NULL) printf("malloc error: steps\n");
	if((steps->end_SNP = (int64_t*) malloc(steps->steps_count * sizeof(int64_t))) == NULL) printf("malloc error: steps\n");
        return steps;
}

void freeSteps(steps_t* steps)
{
	free(steps->start_genetic_position);
	free(steps->end_genetic_position);
	free(steps->start_SNP);
	free(steps->end_SNP);
	free(steps);
}

int64_t countSteps(dbs_t* dbs, filter_pars_t* pars)
{
	return ceil((dbs->genetic_position[dbs->SNP_count - 1] - dbs->genetic_position[0]) / pars->step_size);
}

void setSteps(dbs_t* dbs, steps_t* steps, filter_pars_t* pars)
{
	int64_t i, p = 0;
	int g;
	for(i = 0; i < steps->steps_count; i++)
	{
		if(!i)
		{
			steps->start_genetic_position[0] = dbs->genetic_position[0];
		}
		else
		{
			steps->start_genetic_position[i] = steps->start_genetic_position[i - 1] + pars->step_size;
		}
		steps->end_genetic_position[i] = steps->start_genetic_position[i] + pars->step_size;

		g = 0;
		while(dbs->genetic_position[p] >= steps->start_genetic_position[i] && dbs->genetic_position[p] < steps->end_genetic_position[i])
		{
			if(!g)
			{
				steps->start_SNP[i] = p;
				g = 1;
			}
			steps->end_SNP[i] = p;
			p++;
		}

		if(!g)
		{
			if(!i)
			{
				steps->start_SNP[i] = 0;
				steps->end_SNP[i] = 0;
			}
			else
			{
				steps->start_SNP[i] = steps->end_SNP[i - 1] + 1;
				steps->end_SNP[i] = steps->end_SNP[i - 1];
			}
		}
	}
}

void printSteps(steps_t* steps, dbs_t* dbs)
{
        int64_t i;
        for(i = 0; i < steps->steps_count; i++)
        {
		printf("%f -> %f\t%" PRIu64 " (%f) -> %" PRIu64 " (%f)\n", steps->start_genetic_position[i], steps->end_genetic_position[i], steps->start_SNP[i], dbs->genetic_position[steps->start_SNP[i]], steps->end_SNP[i], dbs->genetic_position[steps->end_SNP[i]]);
        }
}
