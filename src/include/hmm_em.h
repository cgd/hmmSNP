/*
   File Name: hmm_em.h
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
      Algorithm designed by Jin Szatkeiwicz under the supervision of 
      Gary Churchill; and a matlab version of the algorithm is available from 
      Jin.
   Date developed: July-October 2006
   Purpose: struct definitions and function prototype for em algorithm
   Comments:
 
   
   
   
  Copyright (c) 2010 The Jackson Laboratory
 
  This software was developed by Gary Churchill's Lab at The Jackson
  Laboratory (see http://research.jax.org/faculty/churchill).
 
  This is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  This software is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this software. If not, see <http://www.gnu.org/licenses/>.
 */



#ifndef HMMEM_H
#define HMMEM_H
#include "dataIO.h"


#define ELIKE 1


/* define a structure for E-M input parameters */
typedef struct {
   snp_t *snps;
   int num_snps;
   int num_strains;
   
   /* user parameters */
   int max_haplotypes;
   int fix_emission;
   int miss_option;
   int max_iterations;
   int stop_option;
   int partial_output_interval;
    int keep_partial;
   
   int **prune_status_matrix;
   
   double tolerance;
   double ***input_transition_matrix;
   double ***input_emission_matrix;
   
   double lambda_ratio;
   
   double pseudo_option;
   
   char *file_prefix;
   

} em_params_t;


/* structure foro E-M output values */
typedef struct {
   double    aic;
   double    bic[2];
   double ***last_emission_matrix;
   double ***last_transition_matrix;
   double  **last_init_matrix;
   double    last_log_liklihood;
   
   double    elapsed_time;

   int k;
   int num_states;
   int last_iteration;
   int miss;
   int err; /*optional error code */
   
} em_out_t;


em_out_t *hmm_em(em_params_t *input_params, char log_like_out[]);

void freeEM_Out(em_out_t *EM_Out, int num_strains);

double loglike(snp_t *psnpDataArray, int strain, int num_SNPs, 
               double ***emission_matrix, int num_haplotypes, double **pred, 
               int miss_option);
               
              
#endif
