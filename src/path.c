/*
   File Name: path.c
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
      Algorithm designed by Jin Szatkeiwicz under the supervision of 
      Gary Churchill; and a matlab version of the algorithm is available from 
      Jin.
   Date developed: July-October 2006
   Purpose: contains functions to compute the haplotype path
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
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <sys/param.h>



#include "hmm_em.h"
#include "constants.h"
#include "path.h"
#include "filter.h"
#include "hmmUtil.h"
#include "hmmMatrix.h"



/*******************************************************************************
 * maxSmthPath - computes a path by maximizing the smoothness
 *
 * INPUTS: const em_params_t *em_params - input params for EM
 *         const em_out_t *em_out - EM output
 *         const char *prefix - output file prefix
 *         int smooth_out - output final smoothness for each sequence?
 *         int ***strain_haplotype_paths - pointer to 2D integer array to contain path
 *         double ***prob_matrix - pointer to 2D double array to hold 
 *            posterior probability
 *
 * RETURNS: void,  path and probability returne by reference
 *
 * ASSUMES:
 *
 * EFFECTS: allocates memory for path and posterior probability. modifies 
 *    passed by reference pointer to point to the newly allocated memory
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void maxSmthPath(const em_params_t *em_params, const em_out_t *em_out,
                 const char *prefix, int smooth_out,
                 int **strain_haplotype_paths, double **prob_matrix)
{

   FILE *smth_out_fs = NULL;
   char path[MAXPATHLEN + 1];
   char line_buffer[LINE_BUFFER_SIZE];
   
   int max;
 
   
   /* number of rows for Pred,Filt,Smth */
   int num_rows =  em_params->num_snps  + 2;
   
   double **pred_matrix = allocatePred(num_rows, em_params->max_haplotypes);
   double **filt_matrix = allocatePred(num_rows, em_params->max_haplotypes);
   double **smth_matrix = allocatePred(num_rows, em_params->max_haplotypes);
   
   
   for (int i = 0; i < em_params->num_strains; i++)
   {
       filt(pred_matrix, filt_matrix, em_params->snps, i, 
              em_out->last_emission_matrix, em_out->last_transition_matrix, 
              em_out->last_init_matrix, em_params->num_snps, 
              em_params->max_haplotypes, em_params->miss_option);
             
       smooth(em_params->num_snps, em_out->last_transition_matrix, smth_matrix, 
              filt_matrix, pred_matrix, em_params->max_haplotypes, 
              em_params->miss_option);

       /* user may want to save all the smoothness values for each sequence */              
       if (smooth_out)
       {
          /* open file and write out Smth for sequence i */
          strncpy(path, prefix, MAXPATHLEN);
          
          /* XXX - this can end up truncating our suffix
                   we should truncate the prefix if there isn't enough room 
                   in the buffer. non-critical. */
          snprintf(path+strlen(path), MAXPATHLEN - strlen(path), 
                   "smth_%d.csv", i);
          path[MAXPATHLEN-1] = '\0';
  
          if ( (smth_out_fs = fopen(path, "w")) == NULL)
          {
             printf("Unable to open file:\n\t%s\n\t%s\n", path, 
                    strerror(errno));
          }
          else
          {
             line_buffer[LINE_BUFFER_SIZE - 1] = '\0';
             
             /* make a header for the file */
             snprintf(line_buffer, LINE_BUFFER_SIZE - 1, "SNP#,");
             for (int j = 0; j < em_params->max_haplotypes; j++)
             {
                snprintf(line_buffer+strlen(line_buffer), 
                         LINE_BUFFER_SIZE - strlen(line_buffer) - 1,
                         "Hap_%d,", j);
                         
             }
             /* replace trailing ',' with '\n' */
             line_buffer[strlen(line_buffer)-1] = '\n'; 
             /* write string to file */
             fputs(line_buffer, smth_out_fs);
             
             for (int j = 1; j <= em_params->num_snps; j++)
             {
                
                sprintf(line_buffer, "%d,", j);
                for (int k = 0; k < em_params->max_haplotypes; k++)
                {
                   snprintf(line_buffer+strlen(line_buffer), 
                            LINE_BUFFER_SIZE - strlen(line_buffer) - 1, 
                            "%f,", smth_matrix[j][k]);
                }
                /* replace trailing ',' with '\n' */
                line_buffer[strlen(line_buffer)-1] = '\n';
                /* write string to file */
                fputs(line_buffer, smth_out_fs);

                
             }
             fclose(smth_out_fs);
          }
          
          
          
       }
       

       /* construct the haplotype path by maximizing the smoothness */
       for (int j = 1; j <= em_params->num_snps; j++)
       {
          max = 0;
          for (int k = 1; k < em_params->max_haplotypes; k++)
          {
             if (smth_matrix[j][k] > smth_matrix[j][max])
             {
                max = k;
             }
          }
          strain_haplotype_paths[j - 1][i] = max;
          prob_matrix[j - 1][i] = smth_matrix[j][max];          
       } 
       
       
       
   }
   
  freePred(pred_matrix, num_rows);
  freePred(filt_matrix, num_rows);
  freePred(smth_matrix, num_rows);

}
/*******************************************************************************
 * viterbiPath - compute path using viterbi algorithm
 *
 * INPUTS: const em_params_t *em_params - input params for EM
 *         const em_out_t *em_out - output from EM
 *         int ***strain_haplotype_paths - pointer to 2D array to store computed path
 *         double ***log_prob - pointer to 2D array to store
 *           log probability
 *         double **prob_matrix
 *
 * RETURNS: returns void,  
 *
 * ASSUMES:
 *
 * EFFECTS:modifies *strain_haplotype_paths and *log_prob to point to newly 
 *  allocated memory containing path ahd log probability
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void viterbiPath(const em_params_t *em_params, const em_out_t *em_out, 
                 int **strain_haplotype_paths, double **log_prob, 
                 double **prob_matrix)
{
   double stationary_array[em_params->max_haplotypes];
   double **delta; 
   double ***tmp_matrix;
   
   int max;
   
   /* number of rows for Pred,Filt,Smth */
   int num_rows =  em_params->num_snps  + 2;
   
   double **pred_matrix = allocatePred(num_rows, em_params->max_haplotypes);
   double **filt_matrix = allocatePred(num_rows, em_params->max_haplotypes);
   double **smth_matrix = allocatePred(num_rows, em_params->max_haplotypes);
   
   
 
   delta = allocate2Dd(em_params->num_snps, em_params->max_haplotypes);
   
   tmp_matrix = allocateLambda(em_params->num_snps, em_params->max_haplotypes);

   for (int i = 0; i < em_params->max_haplotypes; i++)
   {
      stationary_array[i] = em_out->last_transition_matrix[0][0][i+1];
   }


   for (int iSequence = 0; iSequence < em_params->num_strains; iSequence++)
   {
   
      filt(pred_matrix, filt_matrix, em_params->snps, iSequence, 
             em_out->last_emission_matrix, em_out->last_transition_matrix, 
             em_out->last_init_matrix, em_params->num_snps, 
             em_params->max_haplotypes, em_params->miss_option);
             
      smooth(em_params->num_snps, em_out->last_transition_matrix, smth_matrix, 
             filt_matrix, pred_matrix, em_params->max_haplotypes, 
             em_params->miss_option);
      
      clearLambda(tmp_matrix, em_params->max_haplotypes, em_params->num_snps);
      
   
      for (int i = 0; i < em_params->num_snps; i++)
      {
         for (int j = 0; j < em_params->max_haplotypes; j++)
         {
            delta[i][j] = 0;
         }
      }
   

   
      /* initialize */
      for (int i = 0; i < em_params->max_haplotypes; i++)
      {
         if (em_params->miss_option == MO_RANDOM 
             && em_params->snps[0].iEncodedSNP[iSequence] == 2)
         {
            delta[0][i] = log(stationary_array[i]);
         }
         else
         {
            delta[0][i] = log(stationary_array[i]) 
                         + log(em_out->last_emission_matrix[1][i][em_params->snps[0].iEncodedSNP[iSequence]]);             
         }
      }

      
      /* iterate */
      for (int t = 1; t < em_params->num_snps; t++)
      {


         for (int i = 0; i < em_params->max_haplotypes; i++)
         {
            for (int j = 0; j < em_params->max_haplotypes; j++)
            {
               tmp_matrix[t][i][j] = delta[t - 1][i] 
                                    + log(em_out->last_transition_matrix[t][i][j]);
            }
         }


         
         for (int j = 0; j < em_params->max_haplotypes; j++)
         {
            max = 0;
            for (int k = 1; k < em_params->max_haplotypes; k++)
            {
               if (tmp_matrix[t][k][j] > tmp_matrix[t][max][j])
               {
                  max = k;
               }  
            }
            if (em_params->miss_option == MO_RANDOM 
                && em_params->snps[t].iEncodedSNP[iSequence] == 2)
            {
               delta[t][j] = tmp_matrix[t][max][j];
            }
            else
            {
               delta[t][j] = tmp_matrix[t][max][j] 
                                  + log(em_out->last_emission_matrix[t+1][j][em_params->snps[t].iEncodedSNP[iSequence]]);
            }
 
         }
         
      }
      

      
      /* find last haplotype */
      max = 0;
      for (int i = 1; i < em_params->max_haplotypes; i++)
      {
         if (delta[em_params->num_snps - 1][i] 
             > delta[em_params->num_snps - 1][max])
         {
            max = i;
         } 
      }

      log_prob[em_params->num_snps - 1][iSequence] = 
                             delta[em_params->num_snps-1][max];
      strain_haplotype_paths[em_params->num_snps - 1][iSequence] = max;
      prob_matrix[em_params->num_snps - 1][iSequence] = smth_matrix[em_params->num_snps][max];  
      /* traceback */
      for (int i = em_params->num_snps - 1; i >0; i--)
      {
         max = 0;
         for (int j = 1; j < em_params->max_haplotypes; j++)
         {
            if (tmp_matrix[i][j][strain_haplotype_paths[i][iSequence]] 
                > tmp_matrix[i][max][strain_haplotype_paths[i][iSequence]])
            {
               max = j;
            }
         }
         strain_haplotype_paths[i - 1][iSequence] = max;
         prob_matrix[i - 1][iSequence] = smth_matrix[i][max];  
         log_prob[i - 1][iSequence] = delta[i-1][max] 
                                             - em_out->last_log_liklihood;
      }
      
      
   } /* end sequence loop */
 
   freePred(pred_matrix, num_rows);
   freePred(filt_matrix, num_rows);
   freePred(smth_matrix, num_rows);
}


