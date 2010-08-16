/*
   File Name: prune.c
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
   Date developed:
   Purpose: prune states with low probabilities before running hmmSNP to 
            convergence
   Overview:
   Comments:
      Pruning has been disabled and has not been maintained. This code may 
      not work as is.
   
   
   
   
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
#include <sys/param.h>
 
#include "constants.h"
#include "dataIO.h"
#include "graph.h"
#include "hmm_em.h"
#include "hmmMatrix.h"
#include "hmmUtil.h"
#include "prune.h"

 
 
 
/*******************************************************************************
 * runPruningIterations - run hmm models for pruning purposes, remove (prune) 
 *    states after each model. repeat until no new states are pruned.
 *
 * INPUTS:  em_params - input parameters for hmm_em
 *          snps - array of snps
 *          emission_types - 2 or 3 (number of possible emission types)
 *          prune_rule - which pruning rule we are using
 *          graph_flag - true for graphviz output, false for no output
 *          pruning_option - parameter for pruning rule
 *          summary_file_fs - file stream for summary file output
 *          
 *
 * RETURNS: number of pruning models that were run
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
int runPruningIterations(em_params_t *em_params,  
                         snp_t *snps, int emission_types, 
                         int prune_rule, int graph_flag, double pruning_option,
                         FILE *summary_file_fs)
{

   em_out_t   *em_output;
   double **marg_probability_matrix;
   char path[MAXPATHLEN];
   int model_number = 0;
   int num_pruned;


   snprintf(path, MAXPATHLEN-1, "%sloglke_0.txt", em_params->file_prefix);

   em_output = hmm_em(em_params, path);


   printf("   HMM EM Pruning Model 0 ");
         
   /* check for convergence */
   if (em_output->err == ELIKE)
   {
      printf("    did not converge: loglkelihood decreased at iteration %d\n", 
             em_output->last_iteration);
   }
   else
   {
      printf("    converged at %d\n", em_output->last_iteration);
   }
      
   /* print timing information */
   printf("    E-M Elapsed Time (in iterations): %.2f(s)\n", 
          em_output->elapsed_time);
   printf("    Avg Iteration: %.3f(s)\n", 
          em_output->elapsed_time / em_output->last_iteration);
       

   fprintf(summary_file_fs, SUMMARY_FORMAT, model_number, em_params->tolerance, 
           em_output->last_iteration, em_output->elapsed_time, 
           em_output->num_states, em_output->k, em_output->last_log_liklihood, 
           em_output->aic, em_output->bic[0], em_output->bic[1]);
   fflush(summary_file_fs);

   marg_probability_matrix = allocate2Dd(em_params->num_snps, em_params->max_haplotypes);

   if (marg_probability_matrix == NULL)
   {
      fprintf(stderr, 
           "    Unable to allocate memory for Marginal Probability Matrix\n");
      exit(EHMM_NOMEM);
   }  

   computeMargProb(marg_probability_matrix, em_output->last_transition_matrix, 
                   em_params->max_haplotypes, em_params->num_snps);
                     
   /* write results out to a file */
      
   /* setup filename, if pruning include model number */
   snprintf(path, MAXPATHLEN-1, "%slamboutp_0.csv", em_params->file_prefix);

                
   /* write out the lambda & outp info */
   writeLambdaOutp(path, em_output->last_transition_matrix, 
                   em_output->last_emission_matrix, marg_probability_matrix, 
                   em_params->max_haplotypes, em_params->num_snps, 
                   em_output->miss);
                
   /* write graphviz file */   
   if (graph_flag == TRUE)
   {
      /* setup filename, if pruning include model number */
      snprintf(path, MAXPATHLEN-1, "%sgraph_0.dot", em_params->file_prefix);

      writeGraph(em_output->last_transition_matrix, em_output->last_emission_matrix, 
                 marg_probability_matrix, em_params->num_snps, em_params->max_haplotypes, 
                 snps, em_output->miss, path, 100, 5);
   }
       
  
   /* model zero is the first hmm_em model we ran, now we keep going 
         until we can't prune anymore */            
   model_number = 1;
   num_pruned = prune(em_params->prune_status_matrix, marg_probability_matrix, 
                      snps, em_params->num_snps, 
                      em_params->max_haplotypes, pruning_option);
                 
   /* keep iterating until no more states are pruned */
   while (num_pruned > 0)
   {

   
      /* use output parameters for new input, both point to the same memory */
      /* 
         em_params->input_transition_matrix == em_output->last_transition_matrix
         em_params->input_emission_matrix == em_output->last_emission_matrix 
       */
      
      /* use the prune status matrix to remove pruned states from 
         lambda and outp */
      pruneLambda(em_params->prune_status_matrix, 
                  em_params->input_transition_matrix, em_params->num_snps, 
                  em_params->max_haplotypes, TRUE);
                  
      pruneOutp(em_params->prune_status_matrix, em_params->input_emission_matrix, 
                em_params->num_snps, em_params->max_haplotypes,
                emission_types);
      
      
      /* free em_output (this doesn't fee lambda and outp, but we still 
         have a pointer to them since we are using them as input */
      freeEM_Out(em_output, getNumStrains());
      
      /* setup the path to the next loglke file */
      snprintf(path, MAXPATHLEN-1, "%sloglke_%d.txt", em_params->file_prefix,
               model_number);
      
      /* run hmm_em again, now that we pruned some states */
      em_output = hmm_em(em_params, path);
      

      /* check for convergence */
      if (em_output->err == ELIKE)
      {
         printf("  HMM EM Pruning Model %d did not converge: " 
                "loglkelihood decreased at iteration %d\n", 
                model_number, em_output->last_iteration);
      }
      else
      {
         printf("  HMM EM Pruning Model %d converged at %d\n", 
                model_number, em_output->last_iteration);
      }
            
      /* print timing information */
      printf("    Elapsed Time (in iterations): %.2f(s)\n", 
             em_output->elapsed_time);
      printf("    Avg Iteration: %.3f(s)\n", 
             em_output->elapsed_time / em_output->last_iteration);
        
      /* compute the new marginal state probability */
      computeMargProb(marg_probability_matrix, em_output->last_transition_matrix, 
                      em_params->max_haplotypes, em_params->num_snps); 
                      
      /* write lambda/outp */
      
      /* first get the path ready for the lambda/outpt output file */
      snprintf(path, MAXPATHLEN-1, "%slamboutp_%d.csv", em_params->file_prefix, 
               model_number);
      
      /* and write to a file */
      writeLambdaOutp(path, em_output->last_transition_matrix, 
                      em_output->last_emission_matrix, marg_probability_matrix, 
                      em_params->max_haplotypes, em_params->num_snps, 
                      em_output->miss);
      
      
      /* write the new graphviz file if necessary */
      if (graph_flag == TRUE)
      {
         snprintf(path, MAXPATHLEN-1, "%sgraph_%d.dot", em_params->file_prefix, 
                  model_number);
           
         writeGraph(em_output->last_transition_matrix, 
                    em_output->last_emission_matrix, 
                    marg_probability_matrix, em_params->num_snps, 
                    em_params->max_haplotypes, snps, 
                    emission_types, path, 100, 5);
      }
      /* print model information to the summary file */
      fprintf(summary_file_fs, SUMMARY_FORMAT, model_number, 
              em_params->tolerance, em_output->last_iteration, 
              em_output->elapsed_time, em_output->num_states, em_output->k, 
              em_output->last_log_liklihood, em_output->aic, em_output->bic[0], 
              em_output->bic[1]);
      fflush(summary_file_fs);

      /* look for more states to prune */
      model_number++;
      num_pruned = prune(em_params->prune_status_matrix, marg_probability_matrix, 
                         snps, em_params->num_snps, 
                         em_params->max_haplotypes, pruning_option);
   }
   
   freeEM_Out(em_output, getNumStrains());
   free2Dd(marg_probability_matrix, em_params->num_snps);
   return model_number;

}
 
 
 
/*******************************************************************************
 * prune - remove states with a low probability
 *
 * INPUTS:  int **prune_status_matrix - set to TRUE if state is pruned
 *          double **marg_probability_matrix - marginal state probabilities
 *          snp_t **snps  - SNP data
 *          int num_snps - number of SNPs
 *          int max_haplotypes
 *          double prune_rule - prune threshold
 *
 * RETURNS: number of states pruned
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
int prune(int **prune_status_matrix, double **marg_probability_matrix, snp_t *snps,
          int num_snps, int max_haplotypes, double prune_rule)
{
   int num_pruned = 0;
   
 
   for (int i = 0; i < num_snps; i++)
   {

      for (int j = 0; j < max_haplotypes; j++)
      {
         if (marg_probability_matrix[i][j] <= prune_rule && prune_status_matrix[i][j] == FALSE)
         {
            prune_status_matrix[i][j] = TRUE;
            snps[i].iNumHaplotypes = snps[i].iNumHaplotypes - 1;
            num_pruned++;
         }
      }
   }
 
   return num_pruned;
 
}

/*******************************************************************************
 * pruneLambda - sets lambda values to zero for pruned states
 *
 * INPUTS:  int **prune_status_matrix - TRUE if state is pruned
 *          double ***transition_matrix - lambda matrix
 *          int num_snps
 *          int max_haplotypes
 *          int normalize
 *
 * RETURNS:  void
 *
 * ASSUMES:
 *
 * EFFECTS: sets values of the lambda matrix to zero if the state is pruned
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void pruneLambda(int **prune_status_matrix, double ***transition_matrix, 
                 int num_snps, int max_haplotypes, int normalize)
{

  /* snp1 - snpT-1 H by H*/   
   for (int i = 1; i < num_snps ; i++)
   {  
      for (int j = 0; j < max_haplotypes; j++)
      {
         
         if (prune_status_matrix[i - 1][j] == TRUE)
         {
            for (int k = 0; k < max_haplotypes; k++)
            {
               /* current state trimmed  set row = 0*/
               transition_matrix[i][j][k] = 0.0;
               /* need to set previous column to zero too */
               if (i > 1)
               {
                  transition_matrix[i - 1][k][j] = 0.0;
               }

            } 
            if ( i == 1)
            {
               transition_matrix[0][0][j+1] = 0.0;
            }
         }
      }
      if (normalize == TRUE)
      {
         rowNormalize(transition_matrix[i], max_haplotypes, 
                      max_haplotypes);
         if (i > 1)
         {
            rowNormalize(transition_matrix[i - 1], max_haplotypes, 
                         max_haplotypes);
         }
         else
         {
            rowNormalize(transition_matrix[i - 1], 1, max_haplotypes + 1);
         }
      }
   }
   
   /* snpT H by 1*/
   for (int i = 0; i < max_haplotypes; i++)
   {     
   
      if (prune_status_matrix[num_snps - 1][i] == TRUE)
      {
         transition_matrix[num_snps][i][0] = 0.0;
         
         /* set previous column to zero */
         for (int j = 0; j < max_haplotypes; j++)
         {
            transition_matrix[num_snps - 1][j][i] = 0;
         }
      }
      if (normalize == TRUE)
      {
         rowNormalize(transition_matrix[num_snps - 1], max_haplotypes, 
                      max_haplotypes);
      }
   }

}

/*******************************************************************************
 * pruneOutp - set outp values to zero if state has been pruned
 *
 * INPUTS:  int **prune_status_matrix - TRUE if state is pruned
 *          double ***emission_matrix - outp matrix
 *          int num_snps
 *          int max_haplotypes
 *          int emission_types
 *
 * RETURNS:  void
 *
 * ASSUMES:
 *
 * EFFECTS: sets values of outp to zero if state has been pruned
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void pruneOutp(int **prune_status_matrix, double ***emission_matrix, 
               int num_snps, int max_haplotypes, int emission_types)
{
   for (int i = 1; i <= num_snps; i++)
   {
      for (int j = 0; j < max_haplotypes; j++)
      {

         if (prune_status_matrix[i - 1][j] == TRUE )
         {
            for (int k = 0; k < emission_types; k++)
            {
               /* dOutpMatrix[i][j][k] = dSsyMatrix[i][j][k] / 1; */
               emission_matrix[i][j][k] = 0.0;
            }
         }                 
           
      }
      rowNormalize(emission_matrix[i], max_haplotypes, emission_types);
   }
   
}

/*******************************************************************************
 * rowNormalize - make sure 2d array rows sum to one (normalizes each row)
 *
 * INPUTS:  double **p
 *          int num_rows
 *          int num_columns
 *
 * RETURNS:
 *
 * ASSUMES:
 *
 * EFFECTS:  normalizes each row in the array
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void rowNormalize(double **p, int num_rows, int num_columns)
{
   double row_sum;

   for (int i = 0; i < num_rows; i++)
   {
      row_sum = 0;
      for (int j = 0; j < num_columns; j++)
      {
         row_sum += p[i][j];
      }
      
      if (row_sum == 0)
      {
         row_sum = 1;
      }
      
      for (int j = 0; j < num_columns; j++)
      {
         p[i][j] = p[i][j] / row_sum;
      }
   }
}

/*******************************************************************************
 *  initPruneStatus - allocates and initializes the prune status matrix
 *
 * INPUTS:   int num_snps
 *           int iNumHaplotypes
 *
 * RETURNS: pointer to newly allocated 2d array
 *
 * ASSUMES:
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS: will exit if unable to allocate memory
 *
 * COMMENTS:
 *
 ******************************************************************************/
int **initPruneStatus(int num_snps, int num_haplotypes)
{
   int **prune_status_matrix;
   prune_status_matrix = allocate2Di(num_snps, num_haplotypes);
   
   if (prune_status_matrix == NULL)
   {
      fprintf(stderr, "Unable to allocate PruneStatus\n");
      exit(EHMM_NOMEM);
   }
   
   for (int i = 0; i < num_snps; i++)
   {
      for (int j = 0; j < num_haplotypes; j++)
      {
         prune_status_matrix[i][j] = FALSE; 
      }
   }
   
   return prune_status_matrix;
   
}
