/*
   File Name: hmm_em.c
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz 
      Algorithm designed by Jin Szatkeiwicz under the supervision of 
      Gary Churchill; and a matlab version of the algorithm is available from 
      Jin.
   Date developed: July-October 2006
   Purpose: this file contains the hmm E-M algorithm
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
 
#include <errno.h> 
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/param.h>
#include <time.h>
#include <unistd.h>

#include "CLI.h"
#include "hmmMatrix.h"
#include "constants.h"
#include "dataIO.h"
#include "filter.h"
#include "hmm_em.h"
#include "prune.h"
#include "hmmUtil.h"
#include "temp_file_cleanup.h"

/* internal function prototypes */

               
int countParams(snp_t *snps, int num_snps, int max_haplotypes, 
                int num_emission_types);
                
int getNumStates(snp_t *snps, int num_snps);
               

double logSumPrior(double ***transition_matrix, 
                   double ***transition_pseudoct_matrix, 
                   double ***emission_matrix, int **prune_status_matrix,
                   int num_snps, int max_haplotypes, int num_emission_types, 
                   double pseudo_option);
            
void estep(snp_t *snps, int strain_index, int num_snps, 
           double ***emission_matrix, double ***transition_matrix, 
           int num_haplotypes, double ***sss_matrix, double ***ssy_matrix, 
           double **init_matrix, double **pred_matrix, double **filter_matrix, 
           double **smooth_matrix, int miss_option);

void mstep(double ***transition_matrix, double ***emission_matrix,  
           double ***sss_matrix, double **ss_matrix, double ***ssy_matrix, 
           double **es_matrix, int num_snps, int num_haplotypes, int miss, 
           double ***transition_pseudoct_matrix, int fix_emission, 
           int **prune_status_matrix, double pseudo_option);

/*******************************************************************************
 *  hmm_em - the HMM E-M algorithm.  This is the heart of the program.
 *
 * INPUTS: em_params_ *input_params - pointer to a struct containing all of the 
 *            input parameters
 *         char log_likelihood_out_path - path to loglikelihood output file.
 *
 * RETURNS: pointer to em_out_t struct (see hmm_em.h)
 *
 * ASSUMES: obviously, that it is being passed a pointer to a valid em_params_t
 *   struct
 *
 * EFFECTS:  allocates memory for lambda, outp, etc. If lambda and outp pointers
 *           are not NULL, then it assumes we are providing the initial value.
 *           The memory at these locations will be altered.
 *
 * ERROR CONDITIONS: if unable to allocate memory the program will exit
 *
 * COMMENTS:
 *
 ******************************************************************************/
em_out_t *hmm_em(em_params_t *input_params, char log_likelihood_out_path[])
{


   /* input parameters reference
   typedef struct {
      snp_t *snps;
      int num_snps;
      int num_strains;
   

      int max_haplotypes;
      int iRI_Set;
      int fix_emission;
      int miss_option;
      double tolerance;
      int max_iterations;
      int stop_option;
   
      double ***input_transition_matrix;
      double ***input_emission_matrix;
      int **prune_status_matrix;
    
      double lambda_ratio;
   
   } em_params_t;
   */
   
   double ***transition_matrix;
   double ***transition_pseudoct_matrix;
   double ***emission_matrix;
   double  **init_matrix;
   double  **pred_matrix;
   double  **filter_matrix;
   double  **smooth_matrix; 
   double  **marginal_probability_matrix;
   
   double ***sss_matrix;
   double ***ssy_matrix;
   double **es_matrix;
   double **ss_matrix;
   
   double conv_data;
   double conv_prior;
   double conv_data_prior;
   double log_likelihood;
   double old_likelihood;
   double log_sum_prior;
   double old_log_sum_prior;

   int miss = 2;

   char partial_results_path[MAXPATHLEN + 1];
   char temp_path[MAXPATHLEN + 1];
   
   FILE *log_like_fs = NULL;
   
   clock_t start_time, end_time;
   
   em_out_t *output = malloc(sizeof(em_out_t));
   output->err = 0;
   
   marginal_probability_matrix = NULL;

   if (log_likelihood_out_path != NULL)
   {
      log_like_fs = fopen(log_likelihood_out_path, "w");
      if (log_like_fs == NULL)
      {
         fprintf(stderr, "WARNING: Unable to open loglikelihood output file.\n" 
                 "  path:  %s\n  Error:  %s\n", log_likelihood_out_path, 
                 strerror(errno));
      }
      else
      {
         fprintf(log_like_fs, "iteration\tloglikelihood\tconvergence(data)"
                            "\tlogsumprior\tconvergence(data+prior)\n");
      }
   }
   
   /* setup initial values */
   

   if (input_params->miss_option == MO_EMISSION)
   {
      miss = 3;
   }
   else if (input_params->miss_option == MO_COMPLETE || 
            input_params->miss_option == MO_RANDOM)
   {
      miss = 2;
   }
   output->miss = miss;

   /* setup outp, lambda */                                  


   if (HMM_VERBOSITY > 0) printf("Using passed Lambda Matrix.\n");
   transition_matrix = input_params->input_transition_matrix;
      

   
   /* allocate and initialize lamb_pseudoct */

   transition_pseudoct_matrix = computeLambdaPseudoct(input_params->max_haplotypes, 
                                                 input_params->num_snps, 
                                                 input_params->lambda_ratio,
                                                 input_params->pseudo_option);
   
   /* make sure the pseudo counts are zero for trimmed states */
   pruneLambda(input_params->prune_status_matrix, transition_pseudoct_matrix, 
                 input_params->num_snps, input_params->max_haplotypes, FALSE);
  


   if (HMM_VERBOSITY > 0) printf("Using passed Outp Matrix.\n");
   emission_matrix = input_params->input_emission_matrix;
      
      
      
      
   double stationary_array[input_params->max_haplotypes];
   for(int i = 0; i < input_params->max_haplotypes; i++)
   {
      stationary_array[i] = transition_matrix[0][0][i + 1];
   }

   /* initial prediction */
   if (HMM_VERBOSITY > 0) printf("Allocating Init:  ");
   init_matrix = allocate2Dd(input_params->num_strains, 
                             input_params->max_haplotypes);   
   if (init_matrix != NULL)
   {
      if (HMM_VERBOSITY > 0) printf("Done\n");
   }
   else 
   {
      fprintf(stderr, "Unable to allocate Init\n");
      exit(EHMM_NOMEM);
   }
   
   for (int i = 0; i < input_params->num_strains; i++)
   {
      for (int j = 0; j < input_params->max_haplotypes; j++)
      {
         init_matrix[i][j] = stationary_array[j];
      }
   }



   /* allocate filt and pred TODO - check each malloc call for NULL pointer */
   int first_dim = 
              input_params->num_snps + 2;
    
   if (HMM_VERBOSITY > 0) printf("Allocating Pred Matrix:  ") ; 
   pred_matrix = allocatePred(first_dim, input_params->max_haplotypes);
   if (HMM_VERBOSITY > 0) printf("Done\n");

    
   if (HMM_VERBOSITY > 0) printf("Allocating Filt Matrix:  ");    
   filter_matrix = allocatePred(first_dim, input_params->max_haplotypes);
   if (HMM_VERBOSITY > 0) printf("Done\n");
       

   if (HMM_VERBOSITY > 0) printf("Allocating Smth Matrix:  ");    
   smooth_matrix = allocatePred(first_dim, input_params->max_haplotypes);
   if (HMM_VERBOSITY > 0) printf("Done\n");
         
   /* compute initial likelihood */
   if (HMM_VERBOSITY > 0) printf("Computing Initial likelihood:  ");
   log_likelihood = 0.0;

   for (int i = 0; i < input_params->num_strains; i++)
   {
      if (isTraining(i))
      {
         filt(pred_matrix, filter_matrix, input_params->snps, i,  
              emission_matrix, transition_matrix, init_matrix, 
              input_params->num_snps, input_params->max_haplotypes, 
              input_params->miss_option);
          
         log_likelihood += loglike(input_params->snps, i, 
                             input_params->num_snps, emission_matrix, 
                             input_params->max_haplotypes, pred_matrix, 
                             input_params->miss_option);
      }
   }
   old_likelihood = log_likelihood + 2 * input_params->tolerance;
   if (HMM_VERBOSITY > 0) printf("%f\n", log_likelihood);
   
   if (HMM_VERBOSITY > 0) printf("Computing Initial logsum prior:   ");
   log_sum_prior = logSumPrior(transition_matrix, transition_pseudoct_matrix, 
                              emission_matrix, 
                              input_params->prune_status_matrix, 
                              input_params->num_snps,
                              input_params->max_haplotypes, miss, 
                              input_params->pseudo_option);
     
   old_log_sum_prior = log_sum_prior + 2 * input_params->tolerance;
   if (HMM_VERBOSITY > 0) printf("%f\n", log_sum_prior);

   
   
   /* allocate sss, ssy, es, ss */
   
   
   /* sss has the same dimensions as lambda */
   sss_matrix = allocateLambda(input_params->num_snps,
                               input_params->max_haplotypes);
   
   /* ssy has the same dimensions as outp */
   ssy_matrix = allocateOutp(input_params->num_snps, 
                             input_params->max_haplotypes, miss);
   
   es_matrix = malloc(sizeof(es_matrix[0]) * (input_params->num_snps + 2));
   es_matrix[0] = malloc(sizeof(es_matrix[0][0]));
   for (int i = 1; i <= input_params->num_snps; i++)
   {
      es_matrix[i] = 
           malloc(sizeof(es_matrix[0][0]) * input_params->max_haplotypes);
   }
   es_matrix[input_params->num_snps + 1] = malloc(sizeof(double));

   ss_matrix = malloc(sizeof(ss_matrix[0]) * (input_params->num_snps + 2));
   ss_matrix[0] = malloc(sizeof(ss_matrix[0][0]));
   for (int i = 1; i <= input_params->num_snps; i++)
   {
      ss_matrix[i] = 
           malloc(sizeof(ss_matrix[0][0]) * input_params->max_haplotypes);
   }
   ss_matrix[input_params->num_snps + 1] = malloc(sizeof(ss_matrix[0][0]));




   /* begin EM iterations */
   output->elapsed_time = 0;
   start_time = clock();
   if (HMM_VERBOSITY > 0) printf("Doing EM Iterations...\n");
   for (int iteration = 1; iteration <= input_params->max_iterations; 
        iteration++)
   {
       
      /* initialize sufficient statistics */
      es_matrix[0][0] = 0.0;
      ss_matrix[0][0] = 0.0;
      for (int i = 1; i <= input_params->num_snps; i++)
      {
         for (int j = 0; j < input_params->max_haplotypes; j++)
         {
            es_matrix[i][j] = 0.0;
            ss_matrix[i][j] = 0.0;
         }
      }
      es_matrix[input_params->num_snps + 1][0] = 0.0;
      ss_matrix[input_params->num_snps + 1][0] = 0.0;
      
      clearLambda(sss_matrix, input_params->max_haplotypes,
                  input_params->num_snps);
      clearOutp(ssy_matrix, input_params->max_haplotypes, miss,
                input_params->num_snps);
      
  
      /* for each sequence, run filter and smoother */
      
      for (int i = 0; i < input_params->num_strains; i++)
      {  
         /* if this strain i is to be used for training run the estep on it */
         if (isTraining(i))
         {
            estep(input_params->snps, i, input_params->num_snps, 
                  emission_matrix, transition_matrix, input_params->max_haplotypes, 
                  sss_matrix, ssy_matrix, init_matrix, pred_matrix, filter_matrix, 
                  smooth_matrix, input_params->miss_option);
         }
      } 
 
      
      mstep(transition_matrix, emission_matrix, sss_matrix, ss_matrix, 
            ssy_matrix, es_matrix, input_params->num_snps, 
            input_params->max_haplotypes, miss, 
            transition_pseudoct_matrix, input_params->fix_emission, 
            input_params->prune_status_matrix, input_params->pseudo_option);
         
            
      /* compute log likelihood using updated parameters */
      log_likelihood = 0.0;
      for (int i = 0; i < input_params->num_strains; i++)
      {
         if (isTraining(i))
         {
            filt(pred_matrix, filter_matrix, input_params->snps, i, emission_matrix, 
                   transition_matrix, init_matrix, input_params->num_snps, 
                   input_params->max_haplotypes, input_params->miss_option);
            log_likelihood += loglike(input_params->snps, i, input_params->num_snps, 
                                emission_matrix, input_params->max_haplotypes, 
                                pred_matrix, input_params->miss_option);
         }
      }
           
      log_sum_prior = logSumPrior(transition_matrix, transition_pseudoct_matrix, emission_matrix, 
                              input_params->prune_status_matrix, input_params->num_snps,
                              input_params->max_haplotypes, miss,
                              input_params->pseudo_option);
      
      
      /* check for convergence */
      conv_data = log_likelihood - old_likelihood;
      conv_prior = log_sum_prior - old_log_sum_prior;
      conv_data_prior = (log_likelihood + log_sum_prior) - (old_likelihood + old_log_sum_prior);
      
             
      if (log_like_fs != NULL)
      {
         fprintf(log_like_fs, "%d\t%f\t%f\t%f\t%f\n", iteration, log_likelihood, 
                 conv_data, log_likelihood + log_sum_prior, conv_data_prior);
      }
      
      if (conv_data_prior < 0)
      {
         /* error - likelihood should never decrease */
         output->err = ELIKE;
         output->last_iteration = iteration;
         break;
      }
      
      if (input_params->stop_option == SO_CONV 
          && conv_data_prior < input_params->tolerance)          
      {
         output->last_iteration = iteration;
         break;
      }
      
      old_likelihood = log_likelihood;
      old_log_sum_prior = log_sum_prior;
      output->last_iteration = iteration;

      if (input_params->partial_output_interval > 0 && iteration % input_params->partial_output_interval == 0)
      {
         printf("writing partial output at iteration %d\n", iteration);
         
         if (marginal_probability_matrix == NULL)
         {
            marginal_probability_matrix = allocate2Dd(input_params->num_snps, 
                                            input_params->max_haplotypes);
         }
         
         if (marginal_probability_matrix == NULL)
         {
            fprintf(stderr, 
                 "    Unable to allocate memory for Marginal Probability Matrix\n");
            exit(EHMM_NOMEM);
         }  

         computeMargProb(marginal_probability_matrix, transition_matrix, 
                         input_params->max_haplotypes,
                         input_params->num_snps);
         
         
         if (iteration != input_params->partial_output_interval)
         {
            strcpy(temp_path, partial_results_path);
         }
         
         sprintf(partial_results_path, "%slamboutp_iteration-%d.csv", 
                 input_params->file_prefix, iteration); 
         
         writeLambdaOutp(partial_results_path, transition_matrix, 
                            emission_matrix, marginal_probability_matrix, 
                            input_params->max_haplotypes, 
                            input_params->num_snps, miss);
          
         if (input_params->keep_partial == FALSE && 
             iteration != input_params->partial_output_interval)
         {
            unlink(temp_path);
         }                            
      }
   }
   

   if (input_params->partial_output_interval > 0 &&
       input_params->keep_partial == FALSE &&
       register_temp_file(partial_results_path) != 0)
   {
       printf("Unable to register temp file for automatic cleanup:\n");
       printf("  %s\n", partial_results_path);
       printf("  You may want to delete this file manually.\n");
   }

   
   
   end_time = clock();
   output->elapsed_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC; 
   
   
   output->k = countParams(input_params->snps, input_params->num_snps, 
                            input_params->max_haplotypes, miss);
                            
   output->aic = log_likelihood - output->k;
   output->bic[0] =   log_likelihood - 0.5 * output->k * log(input_params->num_snps);
   output->bic[1] =   log_likelihood - 0.5 * output->k * log(getNumTrainingStrains());  


   /* we return pointers to lamb, outp, init, so we don't free them */
   output->last_transition_matrix = transition_matrix;
   output->last_emission_matrix = emission_matrix;
   output->last_init_matrix = init_matrix;
   output->last_log_liklihood = log_likelihood;
   output->num_states = getNumStates(input_params->snps, 
                                     input_params->num_snps);



   /* Free sss,ssy,pred,filt,smth,es,ss */
   freeLambda(sss_matrix, input_params->num_snps, input_params->max_haplotypes);
   freeOutp(ssy_matrix, input_params->num_snps, input_params->max_haplotypes);
   freePred(pred_matrix, first_dim);
   freePred(filter_matrix, first_dim);
   freePred(smooth_matrix, first_dim);
   for (int i = 0; i < input_params->num_snps + 2; i++)
   {
      free(es_matrix[i]);
   }
   free(es_matrix);

   for (int i = 0; i < input_params->num_snps + 2; i++)
   {
      free(ss_matrix[i]);
   }
   free(ss_matrix);
   
   
   if (marginal_probability_matrix != NULL)
   {
      free2Dd(marginal_probability_matrix, input_params->num_snps);
   }
   
   if (log_like_fs != NULL)
      fclose (log_like_fs);
   
   if (HMM_VERBOSITY > 0) printf("Done EM...\n");
   
   return output;
}



/*******************************************************************************
 * freeEM_Out free an em_out_t type struct and some pointers held by the struct 
 *
 * INPUTS: em_out_t *em_out - a pointer to an em_out_t struct
 *         num_strains - number of strains in input data, needed to fee 
 *             em_out->last_init_matrix
 *
 * RETURNS: void
 *
 * ASSUMES:
 *
 * EFFECTS: frees memory allocated by malloc
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS: the struct also contains pointers to two other multi-dimensional 
 *     arrays that need to be freed, but since these may be allocated elsewhere,
 *     and may be pointed to by other pointers, they will be freed seperately. 
 *
 ******************************************************************************/
void freeEM_Out(em_out_t *em_out, int num_strains)
{
   /* free init, do not free lambda and outp */
   free2Dd(em_out->last_init_matrix, num_strains);
   free (em_out);
}



/*******************************************************************************
 * loglike - compute loglikelihood
 *
 * INPUTS:  snp_t *snps - input data
 *          int strain_index - which sequence we are computing the likelihood for
 *          int num_snps - number of SNPs
 *          double ***emission_matrix - emmissino matrix
 *          int num_haplotypes - number of different happlotypes
 *          double **pred_matrix - pred matrix
 *          int miss_option - MO_RANDOM or MO_EMISSION
 *
 * RETURNS: double loglikelihood value
 *
 * ASSUMES:
 *
 * EFFECTS: none
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
double loglike(snp_t *snps, int strain_index, int num_snps, 
               double ***emission_matrix, int num_haplotypes, double **pred_matrix, 
               int miss_option)
{
   double sum;
   double log_likelihood = 0.0;
   
   for (int i = 0; i < num_snps; i++)
   {
      if ((miss_option == 1 && snps[i].iEncodedSNP[strain_index] < 2) 
          || miss_option == 0 || miss_option == 2)
      {
         sum = 0.0;
         for (int j = 0; j < num_haplotypes; j++)
         {
            sum += emission_matrix[i+1][j][snps[i].iEncodedSNP[strain_index]] * 
                    pred_matrix[i+1][j];
         }
         log_likelihood += log(sum);
      }
      
   }
   return log_likelihood;
}




/* internal functions */



/*******************************************************************************
 *  countParams - count the number of parameters in the hmm model
 *
 * INPUTS:  snp_t *snps
 *          int num_snps
 *          int max_haplotypes, 
 *          int num_emission_types
 *
 * RETURNS: number of parameters
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
int countParams(snp_t *snps, int num_snps, int max_haplotypes, 
                int num_emission_types)
{
   int parameters = 0;
   
   for (int i = 0; i< num_snps; i++)
   { 
       for (int j = 0; j < snps[i].iNumHaplotypes; j++)
       {       
         
            parameters += (num_emission_types - 1);
          
            parameters += (snps[i].iNumHaplotypes - 1);
       }
   }
     
   return parameters;   

}

/*******************************************************************************
 * getNumStates - return the number of states in the model
 *
 * INPUTS:  snp_t *snps
 *          int num_snps
 *
 * RETURNS:  returns number of states in the model
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
int getNumStates(snp_t *snps, int num_snps)
{
   int states = 2; /* start at 2 for begin and end states */
   
   for (int i = 0; i < num_snps; i++)
   {
      states += snps[i].iNumHaplotypes;
   }
     
   return states;
}


/*******************************************************************************
 *  logSumPrior 
 *
 * INPUTS:   double ***transition_matrix  - lambda matrix
 *           double ***transition_pseudoct_matrix - lambda pseudocounts
 *           double ***emission_matrix - emission matrix
 *           int **prune_status_matrix - prune status for each state
 *           int num_snps 
 *           int max_haplotypes
 *           int num_emission_types 
 *           double pseudo_option
 *
 * RETURNS:  returns the sum of the log of all the (lamda/outp values * pseudoptions)
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
double logSumPrior( double ***transition_matrix, 
                    double ***transition_pseudoct_matrix, double ***emission_matrix,
                    int **prune_status_matrix, int num_snps, 
                    int max_haplotypes, int num_emission_types, 
                    double pseudo_option)
 {
 
   double log_sum_prior = 0;
 
   for (int t = 1; t <= num_snps; t++)
   {
      for (int i = 0; i < max_haplotypes; i++)
      {
         if (prune_status_matrix[t - 1][i] == FALSE)      
         {         
            for (int j = 0; j < num_emission_types; j++)
            {
               log_sum_prior += (pseudo_option * get_pseudo_mod())  
                               * log(emission_matrix[t][i][j]);
            }
         }
      }
   }

   for (int t = 1; t < num_snps; t++)
   {
      for (int i = 0; i < max_haplotypes; i++)
      {
         for (int j = 0; j < max_haplotypes; j++)
         {
            if (transition_pseudoct_matrix[t][i][j] != 0)
            {
               
               log_sum_prior += transition_pseudoct_matrix[t][i][j]
                             * log(transition_matrix[t][i][j]);
            }
         }   
      }
   }
   
   for (int i = 0; i < max_haplotypes; i++)
   {
      if (transition_pseudoct_matrix[num_snps][i][0] != 0)
      {
         log_sum_prior += transition_pseudoct_matrix[num_snps][i][0]
                         * log(transition_matrix[num_snps][i][0]);
      }
   }
   
   return log_sum_prior;
}

/*******************************************************************************
 * estep - expectation step . hmm_em has two major steps
 *
 * INPUTS:
 *
 * RETURNS: void, returns sss and ssy by reference
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
void estep(snp_t *snps, int strain_index, int num_snps, 
           double ***emission_matrix, double ***transition_matrix, int num_haplotypes, 
           double ***sss_matrix, double ***ssy_matrix, double **init_matrix, 
           double **pred_matrix, double **filter_matrix, double **smooth_matrix, 
           int miss_option)
{
   filt(pred_matrix, filter_matrix, snps, strain_index, emission_matrix, 
          transition_matrix, init_matrix, num_snps, num_haplotypes, miss_option);
   smooth(num_snps, transition_matrix, smooth_matrix, filter_matrix, pred_matrix, 
          num_haplotypes, miss_option);
          
   for (int i = 0; i < num_haplotypes; i++)
      init_matrix[strain_index][i] = smooth_matrix[1][i];
      
   
   /* add begin and end states */
   filter_matrix[0][0] = 1;
   pred_matrix[0][0] = 1;
   smooth_matrix[0][0] = 1;
   filter_matrix[num_snps + 1][0] = 1;
   pred_matrix[num_snps + 1][0] = 1;
   smooth_matrix[num_snps + 1][0] = 1;

   
   for (int i = 1; i < num_snps; i++)
   {  
      for (int j = 0; j < num_haplotypes; j++)
      {
         for (int k = 0; k < num_haplotypes; k++)
         {
               if (pred_matrix[i+1][k] == 0)
               {
               continue;
               }
         
               sss_matrix[i][j][k] += smooth_matrix[i+1][k] * 
                                        transition_matrix[i][j][k] * 
                                        filter_matrix[i][j] / pred_matrix[i+1][k];
         }
      }
   }
   
   /* End */
   for (int i = 0; i < num_haplotypes; i++)
   {
      if (pred_matrix[num_snps+1][0] == 0)
      {
         continue;
      }
      sss_matrix[num_snps][i][0] += smooth_matrix[num_snps+1][0] * 
                                      transition_matrix[num_snps][i][0] * 
                                      filter_matrix[num_snps][i] / 
                                      pred_matrix[num_snps+1][0];
   }

   /* begin */   
   for (int i = 0; i < num_haplotypes; i++)
   {
      sss_matrix[0][0][i+1] += init_matrix[strain_index][i];
   }
   
   
   /* emission  */
   
   
   for (int i = 0; i < num_snps; i++)
   {              

      if ( (miss_option == MO_RANDOM 
            && snps[i].iEncodedSNP[strain_index] != 2) 
          || miss_option != MO_RANDOM)
      {
         for (int j = 0; j < num_haplotypes; j++)
         { 
            /* access ssy and smooth with i+1 because there is no "begin" in 
               snps, just SNPs */
            /* XXX: we access the SNP_data out of order by looping over the SNP 
               number while we keep the strain number constant, but we access it
               in order elsewhere (iterating over each strain before moving to 
               the next SNP), so we need to take the cache-hit penalty 
               somewhere. */
            ssy_matrix[i+1][j][snps[i].iEncodedSNP[strain_index]] += 
                                                            smooth_matrix[i+1][j];
         }
      }
   }

   
}

/*******************************************************************************
 * mstep - second step of EM  (maximization step)
 *
 * INPUTS:
 *
 * RETURNS: updated value of sss, ss, ssy, es (passed by reference)
 *
 * ASSUMES:
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS: will exit on divide by zero (pseudoption should prevent 
 *   this.
 *
 * COMMENTS:
 *
 ******************************************************************************/
void mstep(double ***transition_matrix, double ***emission_matrix,  
           double ***sss_matrix, double **ss_matrix, double ***ssy_matrix, 
           double **es_matrix, int num_snps, int num_haplotypes, int miss, 
           double ***transition_pseudoct_matrix, int fix_emission, 
           int **prune_status_matrix, double pseudo_option)
{   
     
   /* begin and end states do not emit any SNPs */

   for (int i = 0; i < miss; i++)
   {
      emission_matrix[0][0][i] = 0.0;
      emission_matrix[num_snps+1][0][i] = 0.0;
   }
   
     
  /* add pseudo count; compute row sum at each i */
   for (int i = 1; i <= num_snps; i++)
   { 
     for (int j = 0; j < num_haplotypes; j++)
     {
        ss_matrix[i][j] = 0;
        for (int k = 0; k < miss; k++)
        {
           /* only add pseudo count if state is not pruned */
           if (prune_status_matrix[i - 1][j] == FALSE)
           {
              ssy_matrix[i][j][k] += (pseudo_option * get_pseudo_mod()); 
           }
           ss_matrix[i][j] += ssy_matrix[i][j][k];
        }
               
     } 

   }
              
   /* snp1 - snpT-1 : H by H*/   
   for (int i = 1; i < num_snps ; i++)
   {  
      for (int j = 0; j < num_haplotypes; j++)
      {
         es_matrix[i][j] = 0.0;
         for (int k = 0; k < num_haplotypes; k++)
         {
            /* transition_pseudoct_matrix[i][j][k] == 0 if associated state has 
               been trimmed */
            
            sss_matrix[i][j][k] += transition_pseudoct_matrix[i][j][k];
            es_matrix[i][j] += sss_matrix[i][j][k];
         }
      }
   }
   
   /* snpT : H by 1*/
   for (int i = 0; i < num_haplotypes; i++)
   {
      
      sss_matrix[num_snps][i][0] += transition_pseudoct_matrix[num_snps][i][0];
      es_matrix[num_snps][i] += sss_matrix[num_snps][i][0];
   }


   /* begin 1 by H+1, skip begin[0] */   
   for (int i = 0; i < num_haplotypes; i++)
   {
   
      sss_matrix[0][0][i+1] += transition_pseudoct_matrix[0][0][i+1];
      es_matrix[0][0] += sss_matrix[0][0][i+1];
   }



   if (fix_emission == FALSE)
   {
   
      for (int i = 1; i <=num_snps; i++)
      {
         for (int j = 0; j < num_haplotypes; j++)
         {

            for (int k = 0; k < miss; k++)
            {
               if (ss_matrix[i][j] != 0.0 && prune_status_matrix[i - 1][j] == FALSE)
               {
                  emission_matrix[i][j][k] = ssy_matrix[i][j][k] / ss_matrix[i][j];
               }
               else if (prune_status_matrix[i - 1][j] == TRUE)
               {
                  /* emission_matrix[i][j][k] = ssy_matrix[i][j][k] / 1; */
                  emission_matrix[i][j][k] = 0.0;
               }
               else
               {
                  fprintf(stderr, "divide by zero ss\n");
                  exit(99);
               }
               
            }
            
         }
      }
   }

      

   /* begin 1 by H+1, skip begin[0] */   
   for (int i = 0; i < num_haplotypes; i++)
   {
      transition_matrix[0][0][i+1] = sss_matrix[0][0][i+1] / es_matrix[0][0];
   }

  /* snp1 - snpT-1 H by H*/   
   for (int i = 1; i < num_snps ; i++)
   {  
      for (int j = 0; j < num_haplotypes; j++)
      {
         for (int k = 0; k < num_haplotypes; k++)
         {
            if (es_matrix[i][j] != 0.0 && transition_pseudoct_matrix[i][j][k] != 0)
            {
               transition_matrix[i][j][k] = sss_matrix[i][j][k] / es_matrix[i][j];
            }
            else if (transition_pseudoct_matrix[i][j][k] == 0)
            {
               /* transition_matrix[i][j][k] = sss_matrix[i][j][k] / 1; */
               transition_matrix[i][j][k] = 0.0;
            }
            else
            {
               fprintf(stderr, "divide by zero es\n");
               exit(99);
            }
         }
      }
   }
   
   /* snpT H by 1*/
   for (int i = 0; i < num_haplotypes; i++)
   {     
   
      
      if (es_matrix[num_snps][i] != 0 && transition_pseudoct_matrix[num_snps][i][0] != 0)
      {
         transition_matrix[num_snps][i][0] = sss_matrix[num_snps][i][0] / 
                                         es_matrix[num_snps][i];      
      }
      else if (transition_pseudoct_matrix[num_snps][i][0] == 0)
      {
         transition_matrix[num_snps][i][0] = 0.0;
      }
      else 
      {
         fprintf(stderr, "divide by zero es\n");
         exit(99);
      }

   }


}

