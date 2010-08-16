/*
   File Name: filter.c
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
      Algorithm designed by Jin Szatkeiwicz under the supervision of 
      Gary Churchill; and a matlab version of the algorithm is available from 
      Jin.
   Date developed: July-October 2006
   Purpose: contains filter()  and smooth() functions used by EM
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




#include <stdio.h>
#include <stdlib.h>

#include "constants.h"
#include "dataIO.h"


/*******************************************************************************
 * filt - filt() function used by EM
 *
 * INPUTS: double **pred_matrix 
 *         double **FiltMatrix - filter value - this will be modified 
 *         snp_t **snps - SNP data
 *         int strain_index - which strain we are currently working on
 *         double ***emission_matrix - outp (emission matrix)
 *         double ***transition_matrix - lambda matrix 
 *         double **init_matrix - init
 *         int num_snps- number of SNPs
 *         int num_haplotypes - number of haplotypes
 *         int miss_option - see constants.h for values. specifies if missing as
 *           emission, or missing as random.
 *
 * RETURNS: void,  modifiels FiltMatrix
 *
 * ASSUMES:
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS: will exit if divide by zero occurs. this shouldn't happen
 *      when we use a pseudocount
 *
 * COMMENTS:
 *
 ******************************************************************************/
void filt(double **pred_matrix, double **filt_matrix, snp_t *snps, 
            int strain_index, double ***emission_matrix, double ***transition_matrix, 
            double **init_matrix, int num_snps, int num_haplotypes,
            int miss_option)
{

   double sum;
   
   /* initialize */
   
   for (int i = 0; i < num_haplotypes; i++)
   {
      pred_matrix[1][i] = init_matrix[strain_index][i];
      filt_matrix[1][i] = 0.0;
   }
   
   for (int i = 2; i <= num_snps; i++)
   {
      for (int j = 0; j < num_haplotypes; j++)
      {
         pred_matrix[i][j] = 0.0;
         filt_matrix[i][j] = 0.0;
      }
   }
   

   /* setup initial value */

   if (miss_option == MO_RANDOM && snps[0].iEncodedSNP[strain_index] == 2) 
   {
      for (int j = 0; j < num_haplotypes; j++)
      {
         filt_matrix[1][j] = pred_matrix[1][j];
      }
   }
   else
   {
      sum = 0.0;
      for (int j = 0; j < num_haplotypes; j++)
      {
         sum +=
             emission_matrix[1][j][snps[0].iEncodedSNP[strain_index]] * pred_matrix[1][j];
             
      }
      for (int j = 0; j < num_haplotypes; j++)
      { 
         if (sum == 0.0)
         {
            printf("divide by zero in filter1()\n");

            exit(99);
         }
         filt_matrix[1][j] = emission_matrix[1][j][snps[0].iEncodedSNP[strain_index]] * 
                             pred_matrix[1][j] / sum;
      }
   }

      
   /* filter */
   for (int i = 1; i < num_snps; i++)
   {
      for (int j = 0; j < num_haplotypes; j++)
      {
         for (int k = 0; k < num_haplotypes; k++)
            pred_matrix[i+1][j] += transition_matrix[i][k][j] * filt_matrix[i][k];
      }
         
      if (miss_option == MO_RANDOM && snps[i].iEncodedSNP[strain_index] == 2)
      {
         for (int j = 0; j < num_haplotypes; j++)
            filt_matrix[i + 1][j] = pred_matrix[i + 1][j];
               
      }
      else
      {
         sum = 0.0;
         for (int j = 0; j < num_haplotypes; j++)
         {
            sum += emission_matrix[i+1][j][snps[i].iEncodedSNP[strain_index]] * 
                      pred_matrix[i+1][j];
         }
         for (int j = 0; j < num_haplotypes; j++)
         {
            if (sum == 0.0)
            {
               printf("divide by zero in filter2()\n");
               exit(99);
            }
            filt_matrix[i+1][j] = emission_matrix[i+1][j][snps[i].iEncodedSNP[strain_index]] * 
                                  pred_matrix[i+1][j] / sum;
         }
      }
   }
   


}




/*******************************************************************************
 * smooth - smoothing function used by EM
 *
 * INPUTS: int num_snps - number of SNPs
 *         double ***transition_matrix - lambda (transition) matrix
 *         double **smooth_matrix - matrix to output smoothness
 *         double **filt_matrix - output from filter() function
 *         double **pred_matrix - pred
 *         int num_haplotypes - number of haplotypes
 *         int miss_option - MO_COMPLETE or MO_RANDOM
 *
 * RETURNS: void, modifies matrix at smooth_matrix
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
void smooth(int num_snps, double ***transition_matrix, double **smooth_matrix, 
            double **filt_matrix, double **pred_matrix, int num_haplotypes,
            int miss_option)
{
   /* init */
   for (int i = 1; i <= num_snps; i++)
      for (int j = 0; j < num_haplotypes; j++)
         smooth_matrix[i][j] = 0.0;
          
   for (int i = 0; i < num_haplotypes; i++)
      smooth_matrix[num_snps][i] = filt_matrix[num_snps][i];
      

   /* smooth  */
   for (int i = num_snps-1; i > 0; i--)
   {
      for (int j = 0; j < num_haplotypes; j++)
      {    
         for (int k = 0; k < num_haplotypes; k++)
         {
            if(pred_matrix[i+1][k] == 0)
            {
            continue;
            }
            smooth_matrix[i][j] += transition_matrix[i][j][k] * 
                                  smooth_matrix[i+1][k] / pred_matrix[i+1][k];
         }
         smooth_matrix[i][j] *= filt_matrix[i][j];
      }
   }
   

}

