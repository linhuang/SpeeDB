/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "datastruct.h"
#include "prepDB.h"
#include "filt.h"
#include "dbs.h"
#include "common.h"
#include "precalc.h"

precalc_t* newPrecalc(dbs_t* dbs)
{
        int64_t i;
        precalc_t* precalc = (precalc_t*) malloc(1 * sizeof(precalc_t));
        if(precalc == NULL) printf("malloc error: precalc\n");
        if((precalc->set_size = (int*) malloc(dbs->SNP_count * sizeof(int))) == NULL) printf("malloc error: precalc\n");
        if((precalc->set = (int**) malloc(dbs->SNP_count * sizeof(int*))) == NULL) printf("malloc error: precalc\n");
        for(i = 0; i < dbs->SNP_count; i++)
        {
                if((precalc->set[i] = (int*) malloc(dbs->tuple_count * sizeof(int))) == NULL) printf("malloc error: precalc\n");
        }
        if((precalc->cluster = (int*) malloc(dbs->SNP_count * sizeof(int))) == NULL) printf("malloc error: precalc\n");

        if((precalc->complement_set_size = (int*) malloc(dbs->SNP_count * sizeof(int))) == NULL) printf("malloc error: precalc\n");
        if((precalc->complement_set = (int**) malloc(dbs->SNP_count * sizeof(int*))) == NULL) printf("malloc error: precalc\n");
        for(i = 0; i < dbs->SNP_count; i++)
        {
                if((precalc->complement_set[i] = (int*) malloc(dbs->tuple_count * sizeof(int))) == NULL) printf("malloc error: precalc\n");
        }
        if((precalc->complement_cluster = (int*) malloc(dbs->SNP_count * sizeof(int))) == NULL) printf("malloc error: precalc\n");

        return precalc;
}

void freePrecalc(dbs_t* dbs, precalc_t* precalc)
{
	int64_t i;
	free(precalc->set_size);
	free(precalc->cluster);
	free(precalc->complement_set_size);
	free(precalc->complement_cluster);
	for(i = 0; i < dbs->SNP_count; i++)
	{
		free(precalc->set[i]);
	}
	free(precalc->set);
	for(i = 0; i < dbs->SNP_count; i++)
	{
		free(precalc->complement_set[i]);
	}
	free(precalc->complement_set);
	free(precalc);
}

int similarJaccard(precalc_t* precalc, int i, int j, filter_pars_t* pars)
{
        if((double)(precalc->set_size[i]) < pars->eta_J * (double)(precalc->set_size[j])) return 0;
        if((double)(precalc->set_size[j]) < pars->eta_J * (double)(precalc->set_size[i])) return 0;
        int diff = 0;
        int pointer_i = 0, pointer_j = 0;
        while(1)
        {
                if(pointer_i == precalc->set_size[i])
                {
                        diff += (precalc->set_size[j] - pointer_j);
                        break;
                }
                else if(pointer_j == precalc->set_size[j])
                {
                        diff += (precalc->set_size[i] - pointer_i);
                        break;
                }

                if(precalc->set[i][pointer_i] == precalc->set[j][pointer_j])
                {
                        pointer_i++, pointer_j++;
                }
                else if(precalc->set[i][pointer_i] > precalc->set[j][pointer_j])
                {
                        diff++;
                        pointer_j++;
                }
                else
                {
                        diff++;
                        pointer_i++;
                }
        }
        int A = (precalc->set_size[i] + precalc->set_size[j] - diff) / 2;
        int B = diff + A;
        if((double)(A) / (double)(B) > pars->eta_J) return 1; else return 0;
}

int complementSimilarJaccard(precalc_t* precalc, int i, int j, filter_pars_t* pars)
{
        if((double)(precalc->complement_set_size[i]) < pars->eta_J * (double)(precalc->complement_set_size[j])) return 0;
        if((double)(precalc->complement_set_size[j]) < pars->eta_J * (double)(precalc->complement_set_size[i])) return 0;
        int diff = 0;
        int pointer_i = 0, pointer_j = 0;
        while(1)
        {
                if(pointer_i == precalc->complement_set_size[i])
                {
                        diff += (precalc->complement_set_size[j] - pointer_j);
                        break;
                }
                else if(pointer_j == precalc->complement_set_size[j])
                {
                        diff += (precalc->complement_set_size[i] - pointer_i);
                        break;
                }

                if(precalc->complement_set[i][pointer_i] == precalc->complement_set[j][pointer_j])
                {
                        pointer_i++, pointer_j++;
                }
                else if(precalc->complement_set[i][pointer_i] > precalc->complement_set[j][pointer_j])
                {
                        diff++;
                        pointer_j++;
                }
                else
                {
                        diff++;
                        pointer_i++;
                }
        }
        int A = (precalc->complement_set_size[i] + precalc->complement_set_size[j] - diff) / 2;
        int B = diff + A;
        if((double)(A) / (double)(B) > pars->eta_J) return 1; else return 0;
}

void printPrecalc(dbs_t* dbs, precalc_t* precalc)
{
        int64_t i, j;
        printf("forward\n");
        for(i = 0; i < dbs->SNP_count; i++)
        {
                printf("%d\t", precalc->cluster[i]);
                for(j = 0; j < precalc->set_size[i]; j++)
                {
                        printf("%d ", precalc->set[i][j]);
                }
                printf("\n");
        }

        printf("complement\n");
        for(i = 0; i < dbs->SNP_count; i++)
        {
                printf("%d\t", precalc->complement_cluster[i]);
                for(j = 0; j < precalc->complement_set_size[i]; j++)
                {
                        printf("%d ", precalc->complement_set[i][j]);
                }
                printf("\n");
        }
}

void storePrecalc(dbs_t* dbs, precalc_t* precalc, char* precalc_filename)
{
        FILE* output = fileOpenWB(precalc_filename);

        int64_t i, j;
        for(i = 0; i < dbs->SNP_count; i++)
        {
                fwrite(&precalc->cluster[i], sizeof(int), 1, output);
                fwrite(&precalc->set_size[i], sizeof(int), 1, output);
                for(j = 0; j < precalc->set_size[i]; j++)
                {
                        fwrite(&precalc->set[i][j], sizeof(int), 1, output);
                }
        }

        for(i = 0; i < dbs->SNP_count; i++)
        {
                fwrite(&precalc->complement_cluster[i], sizeof(int), 1, output);
                fwrite(&precalc->complement_set_size[i], sizeof(int), 1, output);
                for(j = 0; j < precalc->complement_set_size[i]; j++)
                {
                        fwrite(&precalc->complement_set[i][j], sizeof(int), 1, output);
                }
        }
        fileClose(output);
}

precalc_t* loadPrecalc(dbs_t* dbs, char* precalc_filename)
{
        FILE* output = fileOpenRB(precalc_filename);
        precalc_t* precalc = newPrecalc(dbs);

        int64_t i, j;
	int ret;
        for(i = 0; i < dbs->SNP_count; i++)
        {
                ret = fread(&precalc->cluster[i], sizeof(int), 1, output);
                ret = fread(&precalc->set_size[i], sizeof(int), 1, output);
                for(j = 0; j < precalc->set_size[i]; j++)
                {
                        ret = fread(&precalc->set[i][j], sizeof(int), 1, output);
                }
        }

        for(i = 0; i < dbs->SNP_count; i++)
        {
                ret = fread(&precalc->complement_cluster[i], sizeof(int), 1, output);
                ret = fread(&precalc->complement_set_size[i], sizeof(int), 1, output);
                for(j = 0; j < precalc->complement_set_size[i]; j++)
                {
                        ret = fread(&precalc->complement_set[i][j], sizeof(int), 1, output);
                }
        }

        fileClose(output);
        return precalc;
}

precalc_t* precalcDBS(dbs_t* dbs, char* precalc_filename, filter_pars_t* pars)
{
        int64_t i, j;
        int64_t max_cluster = -2, max_complement_cluster = -2;
        precalc_t* precalc = newPrecalc(dbs);

        for(i = 0; i < dbs->SNP_count; i++)
        {
                // forward
                precalc->set_size[i] = 0;
                for(j = 0; j < dbs->tuple_count; j++)
                {
                        if(dbs->tuple[j].genotype[i] == 2)
                        {
                                precalc->set[i][precalc->set_size[i]] = j;
                                precalc->set_size[i]++;
                        }
                }

                // complement
                precalc->complement_set_size[i] = 0;
                for(j = 0; j < dbs->tuple_count; j++)
                {
                        if(dbs->tuple[j].genotype[i] == 0)
                        {
                                precalc->complement_set[i][precalc->complement_set_size[i]] = j;
                                precalc->complement_set_size[i]++;
                        }
                }

                // forward
                if(!i)
                {
                        if(precalc->set_size[i] == 0)
                        {
                                precalc->cluster[i] = -1;
                                max_cluster = -1;
                        }
                        else
                        {
                                precalc->cluster[i] = 0;
                                max_cluster = 0;
                        }
                }
                else
                {
                        if(precalc->set_size[i] == 0)
                        {
                                precalc->cluster[i] = -1;
                        }
                        else
                        {
                                for(j = i - 1; j >= 0; j--)
                                {
                                        // check if precalc->set[i] is similar to precalc->set[j]
                                        if(similarJaccard(precalc, i, j, pars))
                                        {
                                                precalc->cluster[i] = precalc->cluster[j];
                                                break;
                                        }
                                }
                                if(j < 0)
                                {
                                        precalc->cluster[i] = max_cluster + 1;
                                        max_cluster++;
                                }
                        }
                }                
		// complement
                if(!i)
                {
                        if(precalc->complement_set_size[i] == 0)
                        {
                                precalc->complement_cluster[i] = -1;
                                max_complement_cluster = -1;
                        }
                        else
                        {
                                precalc->complement_cluster[i] = 0;
                                max_complement_cluster = 0;
                        }
                }
                else
                {
                        if(precalc->complement_set_size[i] == 0)
                        {
                                precalc->complement_cluster[i] = -1;
                        }
                        else
                        {
                                for(j = i - 1; j >= 0; j--)
                                {
                                        // check if precalc->complement_set[i] is similar to precalc->complement_set[j]
                                        if(complementSimilarJaccard(precalc, i, j, pars))
                                        {
                                                precalc->complement_cluster[i] = precalc->complement_cluster[j];
                                                break;
                                        }
                                }
                                if(j < 0)
                                {
                                        precalc->complement_cluster[i] = max_complement_cluster + 1;
                                        max_complement_cluster++;
                                }
                        }
                }
        }
        storePrecalc(dbs, precalc, precalc_filename);

        return precalc;
}

