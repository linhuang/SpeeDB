/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "dbs.h"

dbs_t* newDBS(datasize_t* data_size)
{
        dbs_t* dbs = (dbs_t*) calloc(1, sizeof(dbs_t));
        if(dbs == NULL) printf("malloc error: dbs\n");
        dbs->SNP_count = data_size->SNP_count;
        dbs->haplotype_count = data_size->haplotype_count;
	if(dbs->haplotype_count % 2)
	{
		printf("error: odd haplotype number %" PRIu64 "\n", dbs->haplotype_count);
		exit(1);
	}
        dbs->tuple_count = dbs->haplotype_count / 2;

        int64_t i;
        if((dbs->genetic_position = (double*) malloc(dbs->SNP_count * sizeof(double))) == NULL) printf("malloc error: dbs\n");
        if((dbs->block_SNP = (int*) malloc(dbs->SNP_count * sizeof(int))) == NULL) printf("malloc error: dbs\n");
        if((dbs->step_SNP = (int*) malloc(dbs->SNP_count * sizeof(int))) == NULL) printf("malloc error: dbs\n");
	if((dbs->major_allele = (char*) malloc(dbs->SNP_count * sizeof(char))) == NULL) printf("malloc error: dbs\n");
        if((dbs->minor_allele = (char*) malloc(dbs->SNP_count * sizeof(char))) == NULL) printf("malloc error: dbs\n");
        if((dbs->minor_allele_frequency = (double*) malloc(dbs->SNP_count * sizeof(double))) == NULL) printf("malloc error: dbs\n");
        if((dbs->tuple = (tuple_t*) malloc(dbs->tuple_count * sizeof(tuple_t))) == NULL) printf("malloc error: dbs\n");
        for(i = 0; i < dbs->tuple_count; i++)
        {
                if((dbs->tuple[i].genotype = (int*) calloc(dbs->SNP_count, sizeof(int))) == NULL) printf("malloc error: dbs\n");
        }
	if((dbs->haplotype = (haplotype_t*) malloc(dbs->haplotype_count * sizeof(haplotype_t))) == NULL) printf("malloc error: dbs\n");
        for(i = 0; i < dbs->haplotype_count; i++)
        {
                if((dbs->haplotype[i].allele = (char*) calloc(dbs->SNP_count, sizeof(char))) == NULL) printf("malloc error: dbs\n");
        }
        return dbs;
}

void freeDBS(dbs_t* dbs)
{
	int64_t i;
	free(dbs->genetic_position);
	free(dbs->block_SNP);
	free(dbs->step_SNP);
	free(dbs->major_allele);
	free(dbs->minor_allele);
	free(dbs->minor_allele_frequency);
	for(i = 0; i < dbs->tuple_count; i++)
	{
		free(dbs->tuple[i].genotype);
	}
	free(dbs->tuple);
	for(i = 0; i < dbs->haplotype_count; i++)
	{
		free(dbs->haplotype[i].allele);
	}
	free(dbs->haplotype);
	free(dbs);
}

void storeDBS(char* output_filename, dbs_t* dbs)
{
        int64_t i, j;
        FILE* output = fileOpenWB(output_filename);
        fwrite(&dbs->SNP_count, sizeof(int64_t), 1, output);
        fwrite(&dbs->haplotype_count, sizeof(int64_t), 1, output);
        for(i = 0; i < dbs->SNP_count; i++)
        {
                fwrite(&dbs->genetic_position[i], sizeof(double), 1, output);
		fwrite(&dbs->major_allele[i], sizeof(char), 1, output);
                fwrite(&dbs->minor_allele[i], sizeof(char), 1, output);
                fwrite(&dbs->minor_allele_frequency[i], sizeof(double), 1, output);
                for(j = 0; j < dbs->tuple_count; j++)
                {
                        fwrite(&dbs->tuple[j].genotype[i], sizeof(int), 1, output);
                }
		for(j = 0; j < dbs->haplotype_count; j++)
                {
                        fwrite(&dbs->haplotype[j].allele[i], sizeof(char), 1, output);
                }
        }
	fileClose(output);
}

void printDBS(dbs_t* dbs)
{
        int64_t i, j;
        printf("%" PRIu64 "\t%" PRIu64 "\n", dbs->SNP_count, dbs->haplotype_count);
        for(i = 0; i < dbs->SNP_count; i++)
        {
                if((int)(dbs->genetic_position[i] * 10 + 1e-8) % 10 == 0 && (int)(dbs->genetic_position[i] * 100 + 1e-8) % 10 == 0 && (int)(dbs->genetic_position[i] * 1000 + 1e-8) % 10 == 0)
                {
                        printf("%.0f\t%c\t%c\t%f", dbs->genetic_position[i], dbs->major_allele[i], dbs->minor_allele[i], dbs->minor_allele_frequency[i]);
                }
                else if((int)(dbs->genetic_position[i] * 100 + 1e-8) % 10 == 0 && (int)(dbs->genetic_position[i] * 1000 + 1e-8) % 10 == 0)
                {
                        printf("%.1f\t%c\t%c\t%f", dbs->genetic_position[i], dbs->major_allele[i], dbs->minor_allele[i], dbs->minor_allele_frequency[i]);
                }
                else if((int)(dbs->genetic_position[i] * 1000 + 1e-8) % 10 == 0)
                {
                        printf("%.2f\t%c\t%c\t%f", dbs->genetic_position[i], dbs->major_allele[i], dbs->minor_allele[i], dbs->minor_allele_frequency[i]);
                }
                else
                {
                        printf("%.3f\t%c\t%c\t%f", dbs->genetic_position[i], dbs->major_allele[i], dbs->minor_allele[i], dbs->minor_allele_frequency[i]);
                }
                for(j = 0; j < dbs->tuple_count; j++)
                {
                        printf("\t%d(%c%c)", dbs->tuple[j].genotype[i], dbs->haplotype[j * 2].allele[i], dbs->haplotype[j * 2 + 1].allele[i]);
                }
                printf("\n");
        }

	printf("block SNP\n");
	for(i = 0; i < dbs->SNP_count; i++)
	{
		printf("%f\t%d\t%f\n", dbs->genetic_position[i], dbs->block_SNP[i], dbs->genetic_position[dbs->block_SNP[i]]);
	}

	printf("step SNP\n");
        for(i = 0; i < dbs->SNP_count; i++)
        {
                printf("%f\t%d\t%f\n", dbs->genetic_position[i], dbs->step_SNP[i], dbs->genetic_position[dbs->step_SNP[i]]);
        }
}

dbs_t* loadDBS(char* dbs_filename)
{
        FILE* dbs_file = fileOpenRB(dbs_filename);
        datasize_t* data_size = (datasize_t*) calloc(1, sizeof(datasize_t));
        if(data_size == NULL) printf("malloc error: data_size\n");
	int ret;
        ret = fread(&data_size->SNP_count, sizeof(int64_t), 1, dbs_file);
        ret = fread(&data_size->haplotype_count, sizeof(int64_t), 1, dbs_file);
        printf("SNP count is %" PRIu64 ", haplotype count is %" PRIu64 "\n", data_size->SNP_count, data_size->haplotype_count);
        dbs_t* dbs = newDBS(data_size);

        int64_t i, j;
        for(i = 0; i < dbs->SNP_count; i++)
        {
                ret = fread(&dbs->genetic_position[i], sizeof(double), 1, dbs_file);
		ret = fread(&dbs->major_allele[i], sizeof(char), 1, dbs_file);
                ret = fread(&dbs->minor_allele[i], sizeof(char), 1, dbs_file);
                ret = fread(&dbs->minor_allele_frequency[i], sizeof(double), 1, dbs_file);
                for(j = 0; j < dbs->tuple_count; j++)
                {
                        ret = fread(&dbs->tuple[j].genotype[i], sizeof(int), 1, dbs_file);
                }
		for(j = 0; j < dbs->haplotype_count; j++)
                {
                        ret = fread(&dbs->haplotype[j].allele[i], sizeof(char), 1, dbs_file);
                }
        }

        free(data_size);
        fileClose(dbs_file);

        return dbs;
}

int searchSNP(int SNP_count, double end_position, double *genetic_position)
{
        int lo = 0, hi = SNP_count;
        int mi = (lo + hi) / 2;
        while(1)
        {
                if(end_position == genetic_position[mi])
                {
                        while(1)
                        {
                                if(end_position != genetic_position[mi])
                                {
                                        return mi - 1;
                                }
                                mi++;
                        }
                }
                else if(end_position > genetic_position[mi])
                {
                        lo = mi;
                        mi = (lo + hi) / 2;
                }
                else
                {
                        hi = mi;
                        mi = (lo + hi) / 2;
                }

                if(hi - lo <= 1)
                {
                        if(end_position > genetic_position[hi])
                        {
                                return -1;
                        }
                        else if(end_position < genetic_position[lo])
                        {
                                return lo - 1;
                        }
                        else if(end_position == genetic_position[hi])
                        {
                                return hi;
                        }
                        else
                        {
                                return lo;
                        }
                }
        }
}

void setBlockSNP(dbs_t* dbs, filter_pars_t* pars)
{
        long long i;
        for(i = 0; i < dbs->SNP_count; i++)
        {
                dbs->block_SNP[i] = searchSNP(dbs->SNP_count, dbs->genetic_position[i] + pars->block_size, dbs->genetic_position);
        }
}

void setStepSNP(dbs_t* dbs, filter_pars_t* pars)
{
        long long i;
        for(i = 0; i < dbs->SNP_count; i++)
        {
                dbs->step_SNP[i] = searchSNP(dbs->SNP_count, dbs->genetic_position[i] + pars->step_size, dbs->genetic_position);
        }
}
