/*
   File Name: hmmMatrix.c
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
      Algorithm designed by Jin Szatkeiwicz under the supervision of 
      Gary Churchill; and a matlab version of the algorithm is available from 
      Jin.
   Date developed: July-October 2006
   Purpose: this file contains functions related to various matrixes used 
            in the hmmem function
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


#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <sys/param.h>

#include "CLI.h"
#include "hmmMatrix.h"
#include "constants.h"
#include "hmmUtil.h"
#include "dataIO.h"
#include "rand.h"


/*******************************************************************************
 * randStart - random initialization of starting parameters
 *
 * INPUTS:  em_params_t *emParams - input parameters for hmm_em
 *          int num_starts - number of times to generate the randomized initial
 *             values
 *          char prefix[] - string file output prefix
 * RETURNS: initial lambda and ouput matrixes are returned by pointer in 
 *             emParams
 *          function returns void
 *
 * ASSUMES: assumes emParams has been initialized, except for lambda and outp
 *
 * EFFECTS: allocates memory for lambda and outp
 *
 * ERROR CONDITIONS: will exit with error if unable to malloc memory
 *
 * COMMENTS:
 *
 ******************************************************************************/

void randStart(em_params_t *emParams, int num_starts, char prefix[])
{
   const char *id = "randStart";


   double ***initial_transition_matrix;
   double ***initial_emission_matrix;
   double  **marg_probability_matrix;
   
   double log_likelihood_max = 0;
   int emission_types;
   char path[MAXPATHLEN]; 
   
   char tmp_path[MAXPATHLEN];
   int fds;
   char *tmp_dir;
   
   em_out_t *em_out;

   tmp_dir = getenv("TMPDIR");
   if (tmp_dir != NULL)
   {
      snprintf(tmp_path, MAXPATHLEN, "%s/%slamboutp_tmp.binXXXXXX", tmp_dir, prefix);
   }
   else
   {
       snprintf(tmp_path, MAXPATHLEN, "%s/%slamboutp_tmp.binXXXXXX", DEFAULT_TEMP_DIR, prefix);
   }
   
   printf("TEMP_DIR: %s\n", tmp_path);
   
   fds = mkstemp(tmp_path);
   
   if (fds == -1)
   {
      fprintf(stderr, "%s: Unable to create temp file name.\n", id);
      exit(EHMM);
   }
   
   /* only use mkstemp to avoid race condition between testing for a 
      file's existence and creating the file, so we can just close fds */ 
   close(fds); 

   printf("\n  TEMP FILE NAME: %s\n", tmp_path);
   
   if (emParams->miss_option == MO_EMISSION)
   {
      emission_types = 3;
   }
   else
   {
      emission_types = 2;
   }
   
   emParams->tolerance = 10e-2;
   
   emParams->input_transition_matrix = allocateLambda(emParams->num_snps, 
                                                 emParams->max_haplotypes);
                                        
   if ( emParams->input_transition_matrix == NULL)
   {
      fprintf(stderr, 
              "%s: Unable to allocate memory for Lambda Matrix\n", id);
      exit(EHMM_NOMEM);
   }
   
   emParams->input_emission_matrix = allocateOutp(emParams->num_snps, 
                                             emParams->max_haplotypes,
                                             emission_types);

   
   if (emParams->input_emission_matrix == NULL)
   {
      fprintf(stderr, 
              "%s: Unable to allocate memory for Lambda Matrix\n", id);
      exit(EHMM_NOMEM);
   }
   
   marg_probability_matrix = allocate2Dd(emParams->num_snps, 
                                   emParams->max_haplotypes);
   
   if (marg_probability_matrix == NULL)
   {
      fprintf(stderr, 
              "%s: Unable to allocate memory for Lambda Matrix\n", id);
      exit(EHMM_NOMEM);
   }
   
   initial_transition_matrix = allocateLambda(emParams->num_snps, 
                                            emParams->max_haplotypes);  
   initial_emission_matrix = allocateOutp(emParams->num_snps, 
                                     emParams->max_haplotypes,
                                     emission_types);
                                        
                                        

   
   for (int i = 0; i < num_starts; i++)
   {

                                        
   
      setupInitRand(initial_transition_matrix, initial_emission_matrix, 
                    emParams->max_haplotypes, emission_types, 
                    emParams->num_snps, emParams->num_strains, num_starts, 
                    emParams->lambda_ratio, emParams->fix_emission);
                    

   
      copyLambda(emParams->input_transition_matrix, initial_transition_matrix, 
                 emParams->num_snps, emParams->max_haplotypes);
                 
      copyOutp(emParams->input_emission_matrix, initial_emission_matrix, 
               emParams->num_snps, emParams->max_haplotypes, 
               emission_types);
                   
      computeMargProb(marg_probability_matrix, emParams->input_transition_matrix, 
                      emParams->max_haplotypes, emParams->num_snps);  

 
      em_out = hmm_em(emParams, NULL);

     
      computeMargProb(marg_probability_matrix, em_out->last_transition_matrix, 
                      emParams->max_haplotypes, emParams->num_snps);
    
      if (i == 0 || em_out->last_log_liklihood > log_likelihood_max)
      {
      
         log_likelihood_max = em_out->last_log_liklihood;
                              
         writeLambdaOutpBinary(tmp_path, initial_transition_matrix, 
                           initial_emission_matrix, emParams->max_haplotypes, 
                           emParams->num_snps, emission_types);
      }
     
      
      freeEM_Out(em_out, emParams->num_strains);
    
   }

   

   readLambdaOutpBinary(tmp_path, emParams->input_transition_matrix, 
                        emParams->input_emission_matrix, emParams->max_haplotypes, 
                        emParams->num_snps, emission_types);

   unlink(tmp_path);

   computeMargProb(marg_probability_matrix, emParams->input_transition_matrix, 
                      emParams->max_haplotypes, emParams->num_snps);
                      
   snprintf(path, MAXPATHLEN-1, "%slamboutp_init.csv", prefix);                   
   writeLambdaOutp(path,  emParams->input_transition_matrix, 
                    emParams->input_emission_matrix, marg_probability_matrix, 
                   emParams->max_haplotypes, emParams->num_snps, 
                   emission_types);
                      
   printf("  %s: maxlikelihood = %f\n", id, log_likelihood_max);


   free2Dd(marg_probability_matrix, emParams->num_snps);


              
   freeLambda(initial_transition_matrix, emParams->num_snps, 
              emParams->max_haplotypes);
               
   freeOutp(initial_emission_matrix, emParams->num_snps, 
            emParams->max_haplotypes);         

   
}


/*******************************************************************************
 * setupInitRand - initialized lambda and outp with random values
 * INPUTS:  double ***transition_matrix - pointer to memory allocated for lambda
 *          double ***emission_matrix - pointer to memory allocated for outp
 *          int max_haplotypes - maximum number of haplotypes at each SNP
 *          int emission_types - number of emission types
 *          int num_snps - number of SNPs
 *          int num_strains - number of strains
 *
 * RETURNS: none
 *
 * ASSUMES: memory has been alocated for lambda and oupt
 *
 * EFFECTS: modifies lambda and output
 *
 * ERROR CONDITIONS: none
 *
 * COMMENTS:
 *
 ******************************************************************************/

void setupInitRand(double ***transition_matrix, double ***emission_matrix,
                  int max_haplotypes, int emission_types, 
                  int num_snps, int num_strains, int num_starts, 
                  double lambda_ratio, int fix_emission)
{

   double l[max_haplotypes][max_haplotypes];
   int mm[emission_types];
   
   int num_states;
   double sum;
   int mult = 100;

   
   num_states = 2 + num_snps * max_haplotypes;

   
   assert(lambda_ratio > 0.0 && lambda_ratio < 1.0);
   
   
   get_emission_prior(mm);
   
   
   for (int i = 0; i < max_haplotypes; i++)
   {
      for (int j = 0; j < max_haplotypes; j++)
      {
         if (i == j)
         
         {
            l[i][j] = lambda_ratio;
         }
         else
         {
            l[i][j] = (1 - lambda_ratio) / (max_haplotypes - 1);      
         }
      }
   }

   double tmp = (1 - lambda_ratio) / (max_haplotypes - 1);
   while (tmp * mult < 1)
   {
     mult *= 10;
   }


   /* setup begin state */
   transition_matrix[0][0][0] = 0.0;
   sum = 0.0;
   for (int i = 1; i <= max_haplotypes; i++)
   {
      transition_matrix[0][0][i] = gammaDev(max_haplotypes);
      sum += transition_matrix[0][0][i];
   }
   for (int i = 1; i <= max_haplotypes; i++)
   {
      transition_matrix[0][0][i] = transition_matrix[0][0][i] / sum;
   }
   

   
   /* setup SNPs 1 to T-1 */
   
   if (num_starts == 1)
   {
      for (int i = 1; i <= num_snps - 1; i++)
      {
         for (int j = 0; j < max_haplotypes; j++)
         {
            sum = 0.0;
            for (int k = 0; k < max_haplotypes; k++)
            {
               transition_matrix[i][j][k] = l[j][k];
               sum += l[j][k];
            }

            if (sum == 0.0)
            {
               sum = 1.0;
            }
            for (int k = 0; k < max_haplotypes; k++)
            {
               transition_matrix[i][j][k] = transition_matrix[i][j][k] / sum;
               

            }

         }
      }
   }
   else
   {
      for (int i = 1; i <= num_snps - 1; i++)
      {
      
         for (int j = 0; j < max_haplotypes; j++)
         {
            sum = 0.0;
            for (int k = 0; k < max_haplotypes; k++)
            {
               transition_matrix[i][j][k] = 
                  gammaDev((int)(l[j][k] *   mult));
               sum += transition_matrix[i][j][k];
            }  
            if (sum == 0.0) sum = 1.0;
            for (int k = 0; k < max_haplotypes; k++)
            {
               transition_matrix[i][j][k] = transition_matrix[i][j][k] / sum; 
            }
         }
      }
   }
   
   /* outp SNP 1 to T */
   if (fix_emission == TRUE && max_haplotypes == 2 && emission_types == 2)
   {
      for (int i = 1; i <= num_snps; i++)
      {
         emission_matrix[i][0][0] = 0.99;
         emission_matrix[i][0][1] = 0.01;
         emission_matrix[i][1][0] = 0.01;
         emission_matrix[i][1][1] = 0.99;
      }
   }
   else
   {
      for (int i = 1; i <= num_snps; i++)
      {
         for (int j = 0; j < max_haplotypes; j++)
         {
            sum = 0.0;
            for (int k = 0; k < emission_types; k++)
            {
               emission_matrix[i][j][k] = gammaDev(mm[k]);
               sum += emission_matrix[i][j][k];
            }
         
            if (sum == 0.0) sum = 1.0;
            for (int k = 0; k < emission_types; k++)
            {
               emission_matrix[i][j][k] = emission_matrix[i][j][k] / sum; 

            }
         }
      }
   }
   
   
   /* setup lambda SNP T */
   
   for (int i = 0; i < max_haplotypes; i++)
   {  
      transition_matrix[num_snps][i][0] = 1;
   }
   

   
   
   /* setup being/end for outp */
   for (int i = 0; i < emission_types; i++)
   {
      emission_matrix[0][0][i] = 0.0;
      emission_matrix[num_snps + 1][0][i] = 0.0;
   }
   
}


/*******************************************************************************
 * computeLambdaPseudoct - compute a new lambda matrix that includes a 
 *    pseudocount
 *
 * INPUTS:  
 *          int num_haplotypes
 *          int num_snps
 *          double lambda_ratio - value of the diagonals in lambda 
 *          double pseudo_option - pseudoption value
 *
 * RETURNS: pointer to lambda matrix with pseudocount factored in
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
double ***computeLambdaPseudoct(int max_haplotypes, 
                                int num_snps, double lambda_ratio,
                                double pseudo_option)
{

   double ***transition_pseudo_counts;  

   double l[max_haplotypes][max_haplotypes];
   
  
   assert (lambda_ratio > 0.0 && lambda_ratio < 1.0);

   
   for (int i = 0; i < max_haplotypes; i++)
   {
      /* we initialize lambda to favor the same haplotype  */
      for (int j = 0; j < max_haplotypes; j++)
      {
         if (i == j)
         { 
            l[i][j] = lambda_ratio;
         }
         else
         {
            l[i][j] = (1 - lambda_ratio) / (max_haplotypes - 1);
         }
      }
      
   }
   
   
   transition_pseudo_counts = allocateLambda(num_snps, max_haplotypes);
   
   
   /* setup  b */ 
   transition_pseudo_counts[0][0][0] = 0;
   for (int i = 1; i < max_haplotypes + 1; i++)
   {
      transition_pseudo_counts[0][0][i] = (double)1/max_haplotypes;
   }
   
   /* copy l */
   
   for (int i = 1; i <= num_snps - 1; i++)
   {
      for (int j = 0; j < max_haplotypes; j++)
      {
         for (int k = 0; k < max_haplotypes; k++)
         {
            transition_pseudo_counts[i][j][k] = l[j][k];
         }
      }
   }
   
   /* SNPT */
   
   for (int i = 0; i < max_haplotypes; i++)
   {  
      transition_pseudo_counts[num_snps][i][0] = 1;
   }

   /* Lambda Initial Block Diagonal Matrix is setup */


              
  /* now contribute pseudo count */
   for (int i = 0; i < max_haplotypes + 1; i++)
      transition_pseudo_counts[0][0][i] *= pseudo_option;
   
   for (int i = 1; i < num_snps; i++)
      for (int j = 0; j < max_haplotypes; j++)
         for (int k = 0; k < max_haplotypes; k++)
            transition_pseudo_counts[i][j][k] *= pseudo_option;
            
   for (int i = 0; i < max_haplotypes; i++)
      transition_pseudo_counts[num_snps][i][0] *= pseudo_option;
      
      

   /* printLambda(transition_pseudo_counts, max_haplotypes, num_snps); */



   return transition_pseudo_counts;

}

/*******************************************************************************
 * allocateLambda - allocates memory for lambda
 *
 * INPUTS: int num_snps, int max_haplotypes
 *
 * RETURNS: pointer to lambda
 *
 * ASSUMES:
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:  if malloc fails we could get a segfault
 *
 * COMMENTS:   should check return value of malloc for NULL and print error 
 *    message if malloc fails
 *
 ******************************************************************************/
double ***allocateLambda(int num_snps, int max_haplotypes)
{
   double ***lambda = malloc(sizeof(double**) * (num_snps + 1));
   /* do begin */
   lambda[0] = allocate2Dd(1, max_haplotypes+1);
   
   /* do SNPs */
   for (int i = 1; i < num_snps; i++)
   {
      lambda[i] = allocate2Dd(max_haplotypes, max_haplotypes);
   }
   
   /* do end */
   lambda[num_snps] = allocate2Dd(max_haplotypes, 1);
   
   return lambda;


}

/*******************************************************************************
 * resetlambda - resets lambda to a specific value
 *
 * INPUTS: double ***transition_matrix - pointer to lambda matrix
 *         int num_haplotypes
 *         int num_snps
 *         double v - value to set to
 *
 * RETURNS: void
 *
 * ASSUMES:
 *
 * EFFECTS: modifies lambda
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void resetLambda(double ***transition_matrix, int num_haplotypes, int num_snps, 
                 double v)
{
   for (int i = 0; i < num_haplotypes + 1; i++)
      transition_matrix[0][0][i] = v;
   
   for (int i = 1; i < num_snps; i++)
      for (int j = 0; j < num_haplotypes; j++)
         for (int k = 0; k < num_haplotypes; k++)
            transition_matrix[i][j][k] =v;
            
   for (int i = 0; i < num_haplotypes; i++)
      transition_matrix[num_snps][i][0] = v;
}

/*******************************************************************************
 * clearLambda - sets all elements of lambda to zero
 *
 * INPUTS: double ***transition_matrix - pointer to lambda
 *         int num_haplotypes
 *         int num_snps
 *
 * RETURNS: void
 *
 * ASSUMES:
 *
 * EFFECTS: modifies lambda
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void clearLambda(double ***transition_matrix, int num_haplotypes, int num_snps)
{
   resetLambda(transition_matrix, num_haplotypes, num_snps, 0.0);
}

/*******************************************************************************
 * freeLambda - frees memory used by lambda
 *
 * INPUTS: double ***transition_matrix - pointer to lambda
 *         int num_snps
 *         int num_haplotypes
 *
 * RETURNS: void
 *
 * ASSUMES:
 *
 * EFFECTS: fees memory used by lambda
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void freeLambda(double ***transition_matrix, int num_snps, int num_haplotypes)
{
   /* beginning */
   free2Dd(transition_matrix[0], 1);
   
   for (int i = 1; i < num_snps; i++)
   {
      free2Dd(transition_matrix[i], num_haplotypes);
   }
   
   /*end */
   free2Dd(transition_matrix[num_snps], num_haplotypes);
   
   free(transition_matrix);

}

/*******************************************************************************
 * copyLambda - copies lambda values into new lambda matrix
 *
 * INPUTS: double ***dest_matrix - destination
 *         double ***transition_matrix - source
 *         int num_snps
 *         int num_haplotypes
 *
 * RETURNS: void, copy passed by reference
 *
 * ASSUMES:
 *
 * EFFECTS: modifies destination
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS: uses memcpy when possible to avoid element by element copies
 *
 ******************************************************************************/
void copyLambda(double ***dest_matrix, double ***transition_matrix, int num_snps, 
                int num_haplotypes)
{
   
   /* copy beginning */
   memcpy(dest_matrix[0][0], transition_matrix[0][0], 
          sizeof(double) * (num_haplotypes + 1));
      
   /* copy sub matrices corresponding to SNPs */
   for (int i = 1; i < num_snps; i++)
   {
      for (int j = 0; j < num_haplotypes; j++)
      {
         memcpy(dest_matrix[i][j], transition_matrix[i][j], 
                sizeof(double) * num_haplotypes);
      }
   }
      
   /* copy end - this is out of sequence in memory, so we can't do memcpy 
         -copy one element at a time */   
   for (int i = 0; i < num_haplotypes; i++)
   {
      dest_matrix[num_snps][i][0] = transition_matrix[num_snps][i][0];
   }
   
}


/*******************************************************************************
 * printFullLambda - prints out lambda matrix, including zeros that are not 
 *   stored in our 3-way sustem
 *
 * INPUTS:  double ***lambda - block diagonal lambda
 *          int num_haplotypes - number of haplotypes
 *          int num_snps - numer of snps
 *
 * RETURNS: void
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
void printFullLambda(double ***lambda, int num_haplotypes, int num_snps)
{

   int row_size = num_haplotypes * (num_snps) + 2;
   
   /* print b */
   for (int i = 0; i < num_haplotypes + 1; i++)
      printf("%1.4f  ", lambda[0][0][i]);
      
   for (int i = 0; i < row_size - (num_haplotypes + 1); i++)
      printf("     0  ");
      
   printf("\n");
   
   /* print l */
   for (int i = 0; i < num_snps-1; i++)
   {
      int iOffset = 1 + (i+1)*num_haplotypes;
      
      for (int j = 0; j < num_haplotypes; j++)
      {
         for (int k = 0; k < iOffset; k++)
            printf("     0  ");
            
         for (int k = 0; k < num_haplotypes; k++)
            printf("%1.4f  ", lambda[i+1][j][k]);
            
         
         for (int k = 0; k < row_size - (iOffset + num_haplotypes); k++)
         {
            printf("     0  ");
         }
         printf("\n");

      }
   
   }

   
   for (int i = 0; i < num_haplotypes; i++)
   {
      for (int j = 0; j < row_size - 1; j++)
         printf("     0  ");
      printf("%1.4f \n", lambda[num_snps][i][0]);
      
   }
   
   for (int i = 0; i < row_size; i++)
      printf("     0  ");
   printf("\n");
   
   

}

/*******************************************************************************
 * printLamda - prints lambda, only printing the elements we store
 *
 * INPUTS: double ***lambda - block diagonal lambda
 *         int num_haplotypes - number of haplotypes
 *         int num_snps - number of SNPs
 *
 * RETURNS: void
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
void printLambda(double ***lambda, int num_haplotypes, int num_snps)
{

   
   /* print b */
   for (int i = 0; i < num_haplotypes + 1; i++)
      printf("%1.4f  ", lambda[0][0][i]);
      
      
   printf("\n\n");
   
   /* print l */
   for (int i = 0; i < num_snps-1; i++)
   {
      
      for (int j = 0; j < num_haplotypes; j++)
      {
            
         for (int k = 0; k < num_haplotypes; k++)
            printf("%1.4f  ", lambda[i+1][j][k]);
            
         
         
         printf("\n");
  

      }
             
            printf("\n");
   
   }

   
   for (int i = 0; i < num_haplotypes; i++)
   {
      printf("%1.4f \n", lambda[num_snps][i][0]);
      
   }

   printf("\n");
   
   

}



/*******************************************************************************
 * clearOutp - resets outp matrix to zero
 *
 * INPUTS:  int ***emission_matrix - pointer to outpt matrix
 *          int num_haplotypes - number of haplotypes
 *          int emission_types - 2 or 3 depending on missng as random or emission
 *
 * RETURNS: void
 *
 * ASSUMES: outp has been allocated
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void clearOutp(double ***emission_matrix, int num_haplotypes, int emission_types, 
               int num_snps)
{
   setOutp(emission_matrix, num_haplotypes, emission_types, num_snps, 0.0);
}

/*******************************************************************************
 * setOutp - resets outp matrix 
 *
 * INPUTS:  int ***emission_matrix - pointer to outpt matrix
 *          int num_haplotypes - number of haplotypes
 *          int emission_types - 2 or 3 depending on missng as random or emission
 *          double v - value used to set
 *
 * RETURNS: void
 *
 * ASSUMES: outp has been allocated
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void setOutp(double ***emission_matrix, int num_haplotypes, int emission_types, 
             int num_snps, double v)
{
   for (int i = 0; i < emission_types; i++)
   {
      emission_matrix[0][0][i] = v;
      emission_matrix[num_snps + 1][0][i] = v;
   }
   for (int i = 1; i <= num_snps; i++)
      for (int j = 0; j < num_haplotypes; j++)
         for (int k = 0; k < emission_types; k++)
            emission_matrix[i][j][k] = v;
}


/*******************************************************************************
 * allocateOutp - allocates 3D memory model for outp matrix
 *
 * INPUTS: int num_snps - number of SNPs
 *         int num_haplotypes - number of haplotypes
 *         int emission_types - number of possible emission values (2 or 3)
 *
 * RETURNS: pointer to allocated memory
 *
 * ASSUMES:
 *
 * EFFECTS: allocates memory
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
double ***allocateOutp(int num_snps, int num_haplotypes, int emission_types)
{
    /* allocate */
    double ***emission_matrix = malloc(sizeof(double**) * (num_snps + 2));
    
    emission_matrix[0] = allocate2Dd(1, emission_types);
    emission_matrix[num_snps + 1] = allocate2Dd(1, emission_types);
    
    for (int i = 1; i <= num_snps; i++)
       emission_matrix[i] = allocate2Dd(num_haplotypes, emission_types);

   return emission_matrix;
}

/*******************************************************************************
 * void freeOutp
 *
 * INPUTS: double ***emission_matrix - pointer to outp matrix to deallocate
 *         int num_snps - number of SNPs
 *         int num_haplotypes - max number of haplotypes
 *
 * RETURNS: void
 *
 * ASSUMES: emission_matrix is a valid pointer to appropriately sized data structure
 *
 * EFFECTS: frees all memory associated with emission_matrix
 *
 * ERROR CONDITIONS: bad pointer may cause run time error
 *
 * COMMENTS:
 *
 ******************************************************************************/
void freeOutp(double ***emission_matrix, int num_snps, int num_haplotypes)
{
   free2Dd(emission_matrix[0], 1);
   free2Dd(emission_matrix[num_snps + 1], 1);
   
   for (int i = 1; i <= num_snps; i++)
   {
      free2Dd(emission_matrix[i], num_haplotypes);
   }
}

/*******************************************************************************
 * copyOutp - copies outp value
 *
 * INPUTS: double ***dest_matrix - destination for copy
 *         double ***emission_matrix - source for copy
 *         int num_snps - number of SNPs
 *         int num_haplotypes - number of haplotypes
 *         int emission_types
 *
 * RETURNS: void,  returns copy by reference
 *
 * ASSUMES:
 *
 * EFFECTS: modifies destination
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void copyOutp(double ***dest_matrix, double ***emission_matrix, int num_snps, 
              int num_haplotypes, int emission_types)
{
       
    /* copy begin */
    memcpy(dest_matrix[0][0], emission_matrix[0][0], 
           sizeof(double) * emission_types);
    
    for (int i = 1; i <= num_snps; i++)
       for (int j = 0; j < num_haplotypes; j++)
          memcpy(dest_matrix[i][j], emission_matrix[i][j], 
                 sizeof(double) * emission_types);
    
    /* copy end */
    memcpy(dest_matrix[num_snps + 1][0], emission_matrix[num_snps + 1][0], 
           sizeof(double) * emission_types);

}




/*******************************************************************************
 * allocatePred - allocates memory for pred, filt, smth
 *
 * INPUTS: int first_dim 
 *         int num_haplotypes
 *
 * RETURNS:
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
double **allocatePred(int first_dim, int num_haplotypes)
{
   
   double **pred_matrix = malloc(sizeof(double*) * first_dim);
   pred_matrix[0] = malloc(sizeof(double));
   for (int i = 1; i < first_dim - 1; i++)
   {
      pred_matrix[i] = malloc(sizeof(double) * num_haplotypes);
   }
   pred_matrix[first_dim - 1] = malloc(sizeof(double));

   return pred_matrix;
}  

/*******************************************************************************
 * freePred - frees pred matrix
 *
 * INPUTS:  double **pred_matrix - matrix to free
 *          int row_count number of rows
 *
 * RETURNS:
 *
 * ASSUMES:
 *
 * EFFECTS: frees pred_matrix
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void freePred(double **pred_matrix, int row_count)
{
   for (int i = 0; i < row_count; i++)
   {
      free(pred_matrix[i]);
   }
   free(pred_matrix);
}





/*******************************************************************************
 * computeMargProb - compute the marginal probablilities
 *
 * INPUTS:  double **marg_probability_matrix - marginal state probablilities will be 
 *              stored here
 *          double ***transition_matrix
 *          int num_haplotypes
 *          int num_snps
 *
 * RETURNS:
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
void computeMargProb(double **marg_probability_matrix, double ***transition_matrix, 
                     int num_haplotypes, int num_snps)
{
   
   /* first we need to clear out the old values */
   
   for (int i = 0; i < num_snps; i++)
      for (int j = 0; j < num_haplotypes; j++)
      {
         marg_probability_matrix[i][j] = 0.0;
      }
   
   for (int i = 0; i < num_haplotypes; i++)
   {
      marg_probability_matrix[0][i] = transition_matrix[0][0][i + 1];
   }
   
   for (int i = 1; i < num_snps; i++)
   {
      for (int j = 0; j < num_haplotypes; j++)
      {
         for (int k = 0; k < num_haplotypes; k++)
         {
            marg_probability_matrix[i][k] += 
                       transition_matrix[i][j][k] * marg_probability_matrix[i-1][j];
         }
      }
   }


}

