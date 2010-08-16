/*
   File Name: sort_max_trace.c
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
   Date developed: May-June 2007
   Purpose: state sorting (max trace) algorithms for hmmSNP
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

/* macro to swap two integers */
#define SWAP(a,b)  {      \
                        int _tmp;         \
                        _tmp = a;         \
                        a    = b;         \
                        b = _tmp;         \
                   }


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


#include "dataIO.h"
#include "constants.h"
#include "hmmUtil.h"
#include "hmmMatrix.h"
#include "filter.h"
#include "sort_max_trace.h"


static void relabelRow(double **matrix, int new_labels[], int max_haplotypes);
static void relabelCol(double **matrix, int new_labels[], int max_haplotypes);
static double trace(double **matrix, int iSize);
static void computeNewOrder(int order_array[], int k, int max_haplotypes);
static int factorial(int n);

/*******************************************************************************
 * sort_max_trace - sort using the max trace algorithm
 *
 * INPUTS:  snp_t **snps  - SNP data
 *          double ***transition_matrix
 *          double ***emission_matrix
 *          int num_snps
 *          int max_haplotypes
 *          int miss_option
 *
 * RETURNS: void
 *
 * ASSUMES:
 *
 * EFFECTS: stores sorted lambda in transition_matrix, sorted outp in 
 *    emission_matrix
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void sort_max_trace(snp_t *snps, double ***transition_matrix, 
            double ***emission_matrix, int num_snps, int max_haplotypes, 
            int miss_option)
{
   double *perm_trace_array;
   int **new_labels = allocate2Di(num_snps, max_haplotypes);
   int new_order[max_haplotypes];
   double **tmp1 = allocate2Dd(max_haplotypes, max_haplotypes);
   double **tmp2 = allocate2Dd(max_haplotypes, max_haplotypes);
   int max;
   double new_log_likelihood;
   double log_likelihood;
   double  **pred;
   double  **filt_matrix;
   double  **smth;
   double **init;
   int fact = factorial(max_haplotypes);
  
   perm_trace_array = malloc(sizeof(double) * fact);
   if (perm_trace_array == NULL)
   {
      fprintf(stderr, "Unable to allocate memory in sort_M_T\n");
      exit(EHMM_NOMEM);
   }

   init = allocate2Dd(getNumStrains(), max_haplotypes);   
   if (init == NULL)
   { 
      fprintf(stderr, "Unable to allocate Init\n");
      exit(EHMM_NOMEM); 
   } 

    
   pred = allocatePred(num_snps + 2, max_haplotypes);
   filt_matrix = allocatePred(num_snps + 2, max_haplotypes);
   smth = allocatePred(num_snps + 2, max_haplotypes);
 
   log_likelihood = 0.0;
 
    /* recomputing likelihood for unsorted data */
   for (int i = 0; i < getNumStrains(); i++)
      for (int j = 0; j < max_haplotypes; j++)
         init[i][j] = transition_matrix[0][0][j + 1];
         

   for (int i = 0; i < getNumStrains(); i++)
   {
      filt(pred, filt_matrix, snps, i, emission_matrix, 
             transition_matrix, init, num_snps, 
             max_haplotypes, miss_option);
      log_likelihood += loglike(snps, i, num_snps, 
                          emission_matrix, max_haplotypes, 
                          pred, miss_option);
   }   
   
   
  
 
 
 

   for (int i = 1; i < num_snps; i++)
   {   
      for (int j = 0; j < max_haplotypes; j++)
      {
         new_labels[i][j] = 0;
      }
   }
   
   for (int i = 0; i < max_haplotypes; i++)
   {
      new_labels[0][i] = i;
   }
   
   
   for (int i = 1; i < num_snps; i++)
   {
      for (int j = 0; j < max_haplotypes; j++)
         for (int k = 0; k < max_haplotypes; k++)
         {
            tmp1[j][k] = transition_matrix[i][j][k];
         }
   
      relabelRow(tmp1, new_labels[i-1], max_haplotypes);
      relabelRow(emission_matrix[i], new_labels[i-1], 
                 max_haplotypes);

   
      for (int j = 0; j < fact; j++)
      {
         for (int k = 0; k < max_haplotypes; k++)
            memcpy (tmp2[k], tmp1[k], 
                    sizeof(double)*max_haplotypes); 
      
         computeNewOrder(new_order, j, max_haplotypes);
         relabelCol(tmp2, new_order, max_haplotypes);
         perm_trace_array[j] = trace(tmp2, max_haplotypes);
      }
      
      max = 0;
      for (int j = 1; j < fact; j++)
      {
         if (perm_trace_array[j] > perm_trace_array[max])
         {
            max = j;
         }
      }
      
      /*for (int j = 0; j < max_haplotypes; j++)
      {
         new_labels[i][j] = new_order[j];
      }*/
      
      computeNewOrder(new_labels[i], max, max_haplotypes);

      
      relabelCol(tmp1, new_labels[i], max_haplotypes);
      for (int j = 0; j < max_haplotypes; j++)
      {
         memcpy(transition_matrix[i][j], tmp1[j], 
                sizeof(double) * max_haplotypes);
      }
      
      

   }
   
   /* last SNP */

   relabelRow(transition_matrix[num_snps], new_labels[num_snps-1], 
              max_haplotypes);
   relabelRow(emission_matrix[num_snps], new_labels[num_snps-1], 
              max_haplotypes);

   
   /* compute new loglikelihood */
   new_log_likelihood = 0.0;

   for (int i = 0; i < getNumStrains(); i++)
      for (int j = 0; j < max_haplotypes; j++)
      {
         init[i][j] = transition_matrix[0][0][j + 1];
      }


   for (int i = 0; i < getNumStrains(); i++)
   {
      filt(pred, filt_matrix, snps, i, emission_matrix, 
             transition_matrix, init, num_snps, 
             max_haplotypes, miss_option);
      new_log_likelihood += loglike(snps, i, num_snps, 
                          emission_matrix, max_haplotypes, 
                          pred, miss_option);
   }
   
   printf("Done Sorting.\n  old likelihood: %.15f  new likelihood: %.15f\n", 
           log_likelihood, new_log_likelihood);
   
   if (fabs(new_log_likelihood - log_likelihood) >= 0.000000001)
   {
      printf("likelihood changed - bug in sorting algorithm!\n");
      exit(EHMM_SORT); 
   
   }
  
   freePred(pred, num_snps + 2);
   freePred(filt_matrix, num_snps + 2);
   freePred(smth, num_snps + 2);
   free(perm_trace_array);
   
}


/*******************************************************************************
 * relabelRow - swaps rows of a 2D array around
 *
 * INPUTS:  double **matrix - 2D array
 *          int new_labels[] - array containing new order of rows
 *          int max_haplotypes
 *
 * RETURNS:  void
 *
 * ASSUMES:
 *
 * EFFECTS:  reorders the rows in matrix based on the seriese of row numbers in 
 *    new_labels
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
static void relabelRow(double **matrix, int new_labels[], int max_haplotypes)
{
   double *tmp[max_haplotypes];

   
   for (int i = 0; i < max_haplotypes; i++)
   { 
      tmp[i] = matrix[i];
   }
   
   for (int i = 0; i < max_haplotypes; i++)
   {
      assert(new_labels[i] < max_haplotypes);
      matrix[i] = tmp[new_labels[i]];
   }
   

}

/*******************************************************************************
 * relabelCol - swaps columns of 2D array around
 *
 * INPUTS:  double **matrix - 2D array
 *          int new_labels[] - new column order
 *          int max_haplotypes
 *
 * RETURNS: void
 *
 * ASSUMES:
 *
 * EFFECTS: swaps columns around based on new column order in new_labels
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
static void relabelCol(double **matrix, int new_labels[], int max_haplotypes)
{
   double tmp[max_haplotypes][max_haplotypes];
   
   for (int i = 0; i < max_haplotypes; i++)
   {
      memcpy(tmp[i], matrix[i], max_haplotypes * sizeof(double));
   }
   
   for (int i = 0; i < max_haplotypes; i++)
   {
      for (int j = 0; j < max_haplotypes; j++)
      {
         matrix[i][j] = tmp[i][new_labels[j]];
      }
   }
}

/*******************************************************************************
 * trace - computes trace of 2D matrix
 *
 * INPUTS: double **matrix - 2D matrix
 *         int iSize -size of square 2D matrix
 *
 * RETURNS: trace value
 *
 * ASSUMES:
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
static double trace(double **matrix, int size)
{
   double trace = 0.0;
     
  for (int i = 0; i < size; i++)
  {
      trace += matrix[i][i];
  }
   
   return trace;

}

/*******************************************************************************
 * computeNewOrder - computes a unique permutation of all numbers 0 through
 *     max_haplotypes for each value of k, k <= total # of permutations
 *
 * INPUTS: int order_array[] - new permutation will be stored here
 *         int k - this deteremines which permutation is generated
 *         int max_haplotypes
 *
 * RETURNS: void
 *
 * ASSUMES:
 *
 * EFFECTS:  alters order_array[]
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
static void computeNewOrder(int order_array[], int k, int max_haplotypes)
{
   int fact = 1;
   for (int i = 0; i < max_haplotypes; i++)
   {
      order_array[i] = i;
   }
   for (int i = 2; i <= max_haplotypes; i++)
   {
      fact *= i-1;
      SWAP(order_array[i - ((k / fact) % i) - 1],order_array[i -1]);
   }

}



/*******************************************************************************
 * factorial - iterative factorial function
 *
 * INPUTS: int n - number to compute factorial of
 *
 * RETURNS: n!
 *
 * ASSUMES:
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS: going with the iteratvie implementation instead of recursive to 
 *   reduce function call overhead...
 *
 ******************************************************************************/
static int factorial(int n)
{
   int result = 1;
   for (int i = 2; i <= n; i++)
   {
      result *= i;
   }
   return result;
}
