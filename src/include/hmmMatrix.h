/*
   File Name: hmmMatrix.h
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
      Algorithm designed by Jin Szatkeiwicz under the supervision of 
      Gary Churchill; and a matlab version of the algorithm is available from 
      Jin.
   Date developed: July-October 2006
   Purpose: function prototypes for externally referenced functions defined in 
      matrix.c
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

#ifndef HMM_MATRIX_H
#define HMM_MATRIX_H

/* our block diagonal lambda matrix is stored without the large zero portions of the 
   matrix.  We store it as a number of smaller sub-matrices */


#include "dataIO.h"
#include "hmm_em.h"

void randStart(em_params_t *em_params, int num_starts, char prefix[]);
void setupInitRand(double ***transition_matrix, double ***emission_matrix,
                  int max_haplotypes, int emission_types, int num_snps, 
                  int num_strains, int num_starts, double lambda_ratio,
                  int fix_emission);


/* lambda related functions */
double ***computeLambdaPseudoct(int max_haplotypes, int num_snps, 
                                double lambda_ratio, double pseudo_option);
double ***allocateLambda(int num_snps, int max_haplotypes);
void resetLambda(double ***transition_matrix, int max_haplotypes, int num_snps, 
                 double v);
void clearLambda(double ***lambda, int max_haplotypes, int num_snps);
void freeLambda(double ***lambda, int num_snps, int max_haplotypes);
void copyLambda(double ***dest, double ***lambda, int num_snps, 
                int max_haplotypes);
void printFullLambda(double ***lambda, int max_haplotypes, int num_snps);
void printLambda(double ***lambda, int max_haplotypes, int num_snps);

/* outp related functions */

double ***allocateOutp(int num_snps, int max_haplotypes, int emission_types);
void freeOutp(double ***emission_matrix, int num_snps, int max_haplotypes);
void copyOutp(double ***dest, double ***emission_matrix, int num_snps, 
              int max_haplotypes, int emission_types);
void setOutp(double ***emission_matrix, int max_haplotypes, int emission_types, 
             int num_snps, double v);
void clearOutp(double ***emission_matrix, int max_haplotypes, int emission_types, 
               int num_snps);

/* allocate pred,filt,smth */
double **allocatePred(int first_dim, int max_haplotypes);
void freePred(double **pred_matrix, int first_dim);


void computeMargProb(double **marg_prob_matrix, double ***transition_matrix, 
                     int max_haplotypes, int num_snps);



#endif
