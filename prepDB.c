/*
 * Lin Huang <linhuang@cs.stanford.edu> 2014
 */

#include <stdlib.h>
#include <string.h>
#include "datastruct.h"
#include "prepDB.h"
#include "common.h"
#include "dbs.h"

datasize_t* countSNPs(char* input_filename)
{
        datasize_t* data_size = (datasize_t*) calloc(1, sizeof(datasize_t));
        if(data_size == NULL) printf("malloc error: data_size\n");

	FILE* input = fileOpenR(input_filename);

	int i;
	char c;
	while((c = (char) getc(input)) != EOF)
	{
		for(i = 0; i < TPED_HEADER; i++)
		{
			while(c != ' ' && c != '\t')
			{
				c = (char) getc(input);
			}

			while(c == ' ' || c == '\t')
			{
				c = (char) getc(input);
			}
		}
		if(c != ' ' && c != '\t')
                {
                        if(data_size->SNP_count == 0) data_size->haplotype_count++;
                }

		while(c != '\n')
		{
			c = (char) getc(input);
			if(c != ' ' && c != '\t' && c != '\n')
			{
				if(data_size->SNP_count == 0) data_size->haplotype_count++;
			}
		}
		data_size->SNP_count++;
	}
	fileClose(input);

        return data_size;
}

void setAllele(char c, char* allele, char* allele_alt, int64_t* allele_occ, int64_t* allele_alt_occ)
{
	if(!isMissing(c))
	{
		if(c == (*allele)) (*allele_occ)++;
		else if(c == (*allele_alt)) (*allele_alt_occ)++;
		else if((*allele) == ' ') (*allele) = c, (*allele_occ)++;
		else if((*allele_alt) == ' ') (*allele_alt) = c, (*allele_alt_occ)++;
	}
}

dbs_t* identifyAlleles(char* input_filename, datasize_t* data_size)
{
	dbs_t* dbs = newDBS(data_size);
        FILE* input = fileOpenR(input_filename);

	int64_t i;
	int64_t SNP_index = 0;
	char gen_pos[MAX_GEN_POS_LENGTH];
	char allele, allele_alt;
	int64_t allele_occ, allele_alt_occ;
	char c;
	while((c = (char) getc(input)) != EOF)
        {
		allele = ' ';
		allele_alt = ' ';
		allele_occ = 0;
		allele_alt_occ = 0;
		memset(gen_pos, 0, MAX_GEN_POS_LENGTH * sizeof(char));
		for(i = 0; i < TPED_HEADER; i++)
                {
			while(c != ' ' && c != '\t')
                        {
				if(i == 2) gen_pos[strlen(gen_pos)] = c;
				c = (char) getc(input);
                        }

                        while(c == ' ' || c == '\t')
                        {
                                c = (char) getc(input);
                        }
		}
		dbs->genetic_position[SNP_index] = atof(gen_pos);
		if(SNP_index >= 1 && dbs->genetic_position[SNP_index] < dbs->genetic_position[SNP_index - 1]) printf("warning: the %" PRIu64 "th marker's genetic position (%f) is smaller than the %" PRIu64 "th marker's (%f)\n", SNP_index, dbs->genetic_position[SNP_index], SNP_index - 1, dbs->genetic_position[SNP_index - 1]);
		setAllele(c, &allele, &allele_alt, &allele_occ, &allele_alt_occ);
		
		while(c != '\n')
                {
                        c = (char) getc(input);
                        if(c != ' ' && c != '\t' && c != '\n')
                        {
				setAllele(c, &allele, &allele_alt, &allele_occ, &allele_alt_occ);
			}
		}
		if(allele_occ >= allele_alt_occ)
		{
			dbs->major_allele[SNP_index] = allele;
			dbs->minor_allele[SNP_index] = allele_alt;
			dbs->minor_allele_frequency[SNP_index] = (double)allele_alt_occ / (double)dbs->haplotype_count;
		}
		else
		{
			dbs->major_allele[SNP_index] = allele_alt;
			dbs->minor_allele[SNP_index] = allele;
			dbs->minor_allele_frequency[SNP_index] = (double)allele_occ / (double)dbs->haplotype_count;
		}
		if(dbs->major_allele[SNP_index] == ' ') 
		{
			dbs->major_allele[SNP_index] = 'N';
			printf("warning: cannot find any valid allele in marker %" PRIu64 "\n", SNP_index);
		}
		if(dbs->minor_allele[SNP_index] == ' ') dbs->minor_allele[SNP_index] = 'N';

		SNP_index++;
	}

	fileClose(input);
        return dbs;
}

void computeDBS(char* input_filename, datasize_t* data_size, dbs_t* dbs)
{
        FILE* input = fileOpenR(input_filename);

	int64_t i;
        char c;
	int64_t haplotype_index = 0, SNP_index = 0;
	while((c = (char) getc(input)) != EOF)
        {
                for(i = 0; i < TPED_HEADER; i++)
                {
                        while(c != ' ' && c != '\t')
                        {
                                c = (char) getc(input);
                        }

                        while(c == ' ' || c == '\t')
                        {
                                c = (char) getc(input);
                        }
                }
		dbs->haplotype[haplotype_index].allele[SNP_index] = c;
		haplotype_index++;

                while(c != '\n')
                {
                        c = (char) getc(input);
                        if(c != ' ' && c != '\t' && c != '\n')
                        {
				dbs->haplotype[haplotype_index].allele[SNP_index] = c;
				if(haplotype_index % 2)
				{
					if(isMissing(dbs->haplotype[haplotype_index - 1].allele[SNP_index]) || isMissing(dbs->haplotype[haplotype_index].allele[SNP_index]))
					{
						dbs->tuple[haplotype_index / 2].genotype[SNP_index] = 1;
					}
					else
					{
						if(dbs->haplotype[haplotype_index - 1].allele[SNP_index] == dbs->major_allele[SNP_index]) dbs->tuple[haplotype_index / 2].genotype[SNP_index]++;
						if(dbs->haplotype[haplotype_index].allele[SNP_index] == dbs->major_allele[SNP_index]) dbs->tuple[haplotype_index / 2].genotype[SNP_index]++;
					}
				}
				haplotype_index++;
                        }
                }

                SNP_index++;
		haplotype_index = 0;
        }

	fileClose(input);
}

int genDB(char* input_filename, char* output_filename)
{
	datasize_t* data_size = countSNPs(input_filename);
	printf("SNP count is %" PRIu64 ", haplotype count is %" PRIu64 "\n", data_size->SNP_count, data_size->haplotype_count);
	dbs_t* dbs = identifyAlleles(input_filename, data_size);
	computeDBS(input_filename, data_size, dbs);
	storeDBS(output_filename, dbs);

	free(data_size);
	freeDBS(dbs);

	return 0;
}
