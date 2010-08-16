/*
   File Name: fill.c
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
      Algorithm designed by Jin Szatkeiwicz under the supervision of 
      Gary Churchill; and a matlab version of the algorithm is available from 
      Jin.
   Date developed: July-October 2006
   Purpose: contains functions used to fill-in data
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
#include <ctype.h>
#include <assert.h>
#include <time.h>


#include "constants.h"
#include "hmmUtil.h"
#include "dataIO.h"
#include "fill.h"
#include "path.h"


/*******************************************************************************
 * fill_and_save - fill in missing genotypes and save the result
 *
 * INPUTS: snp_t *snps - array of snp_t structs
 *         em_params_t *em_params - input parameters for hmm_em()
 *         em_out_t *em_out - output of hmm_em()
 *         int path_option - which path algorithm was used
 *         int sort_option - state sorting option (for path output)
 *         double threshold - confidence threshold for filtered output
 *         int num_emission_types - 2 or 3 depending if missing is a valid emission type
 *         int chromosome - which chromosome was processed
 *         int smooth_out - boolean value to instruct maxSmthPath to write out
 *             the smoothness information
 *         double **fill_probability_matrix - pointer to memory to store filling probabilities 
 *         
 *
 * RETURNS:
 *      funciton returns void, filling probabilities returned via parameter
 *
 * ASSUMES:
 *
 * EFFECTS:  modifies contents of fill_probability_matrix, writes filled file and path file
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void fill_and_save(snp_t *snps, em_params_t *em_params, em_out_t *em_out, 
                   int path_option, int sort_option, double threshold,
                   int num_emission_types, int chromosome, int smooth_out,
                   double **fill_probability_matrix)
{
  clock_t tStart;
  clock_t tEnd;
  int **path_matrix;
  double **log_prob_matrix = NULL;

  assert(path_option == PATH_VITERBI || path_option == PATH_MAXSMTH);
  
  path_matrix = allocate2Di(em_params->num_snps, em_params->num_strains);

  tStart = clock();
  
  if (path_option == PATH_VITERBI)
  {
    printf("Calculating Viterbi path...");
    log_prob_matrix = allocate2Dd(em_params->num_snps, em_params->num_strains);

    viterbiPath(em_params, em_out, path_matrix, log_prob_matrix,
               fill_probability_matrix);
  }
  else
  {
    printf("Calculating max smooth path...");
    maxSmthPath(em_params, em_out, em_params->file_prefix, smooth_out, 
                path_matrix, fill_probability_matrix);
  }
  
  tEnd = clock();
  
  printf("Done (%.2fs)\n", ((double)(tEnd - tStart)) / CLOCKS_PER_SEC);
          
  /* write the path to a file */                
  writePath(snps, path_matrix, em_params->num_snps, 
            em_params->num_strains, em_params->file_prefix, 
            path_option, chromosome, sort_option, 
            em_out->last_emission_matrix);
  
  /* if random missing then fill in using path */
  if ( em_params->miss_option > 0)
  {
    fillin(em_params->snps,  path_matrix, em_out->last_emission_matrix, 
           em_params->num_snps, em_params->num_strains, num_emission_types, 
           fill_probability_matrix);
  
    writeFilledData(snps, chromosome, fill_probability_matrix,
                    em_params->file_prefix, path_option, sort_option, 
                    threshold, em_params->num_snps);
  }

  free2Di(path_matrix, em_params->num_snps);
  if (path_option == PATH_VITERBI)
  {
    free2Dd(log_prob_matrix, em_params->num_snps);
  }
}

/*******************************************************************************
 * fillin - fill in encoded data (data encoded in the 0,1,2 trinary format )
 *
 * INPUTS: snp_t *snps - encoded data
 *         int **piMathMatrix - path computed from viterbi or maxsmth algorithms
 *         double ***emission_matrix - outp
 *         int num_snps - number of SNPs
 *         int num_strains - number of strains
 *         int miss - number of possible output values
 *         double **probability_matrix - fill in probabilities
 *
 * RETURNS: a pointer to a 2D array of fill in probabilities
 *
 * ASSUMES:
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:  We fill in missing (2) with 3 and 4 (3 = major allele, 
 *            4 for minor allele).  This will allow us to fill the data twice
 *            without copying it first (because of the 3 and 4s, we can 
 *            "roll back" the fill in, and re fill with a different algorithm).
 *  
 ******************************************************************************/
void fillin(snp_t *snps, int **path_matrix, double ***emission_matrix, 
            int num_snps, int num_strains, int miss, double **probability_matrix)
{
   int fillin;
   int haplotype;

   

   for (int i = 0; i < num_snps; i++)
   {
      for (int j = 0; j < num_strains; j++)
      {

         /* if data is missing, or has been filled in with 3 or 4 (see function
             header comment) */
         if (snps[i].iEncodedSNP[j] >= 2)
         {
            /* find the most likely fillin value for a given haplotype */
            fillin = 0;
            haplotype = path_matrix[i][j];
            for (int k = 1; k < miss; k++)
            {
               if (emission_matrix[i + 1][haplotype][k] > 
                   emission_matrix[i + 1][haplotype][fillin])
               {
                  fillin = k;
               }
            }
            probability_matrix[i][j] *= 
                                emission_matrix[i + 1][path_matrix[i][j]][fillin];
            if (fillin != 2)
            {
               /* filled zeros are written as 3, filled ones are written as 4 */
               snps[i].iEncodedSNP[j] = fillin + 3; 
            }
            else
            {
               snps[i].iEncodedSNP[j] = 5; /* missing as emission */
            }   
         }
         else
         {
            /* if there was already a zero or one we just copy it */
            probability_matrix[i][j] *= 
                emission_matrix[i + 1][path_matrix[i][j]][snps[i].iEncodedSNP[j]];
         }       
      }
   }
   

}

/*******************************************************************************
 * calculate_fill_rule - calculate the "fill in rules", that is which 
 *      nucleotide each haplotype will be filled with
 *
 * INPUTS:   snp_t snp - the SNP we will calculate the fill rule for
 *           double **emission_for_snp - subset of emission matrix corresponding
 *              to this SNP
 *           int haplotypes - number of haplotypes
 *           int emission_types - number of possible emission types
 *           fill_rule[] - character array to hold the fill rules
 *
 * RETURNS:  function returns void, fill returns are stored in fill_rule[] 
 *            (including terminating NULL character)
 *
 * ASSUMES:  fill_rule[] contains space for haplotypes+1 characters
 *
 * EFFECTS:  modifies fill_rule[]
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void calculate_fill_rule(snp_t snp, double **emission_for_snp, 
                         int haplotypes, int emission_types, char fill_rule[])
{

   int fill;
   fill_rule[haplotypes] = '\0';
   for (int i = 0; i < haplotypes; i++)
   {
      fill = 0;
      for (int j = 1; j < emission_types; j++)
      {
         if (emission_for_snp[i][j] > 
             emission_for_snp[i][fill])
         {
            fill  = j;
         }
      }
      if (emission_for_snp[i][fill] == 0.0)
      {
         fill_rule[i] = '-';
      }
      else if (fill < 2)
      {
         fill_rule[i] = tolower(snp.cAlleleArray[fill]);
      }
      else
      {
         /* missing as emission */
         fill_rule[i] = 'n';
      }
   }
}

/*******************************************************************************
 XXX
 * fillSNP - takes a SNP string and fills it in using the fille encoded SNP 
 *  as a template
 *
 * INPUTS: char snp_str[] - string representation of SNP that needs to be filled
 *         snp_t *snp_enc - encoded SNP and allele information
 *         int num_strains - lenght of SNP
 *         
 * RETURNS: void
 *
 * ASSUMES:
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:  a filled zero is encoded as a 3, and a filled one as a 4:
 *         see fillin()
 *  
 ******************************************************************************/
void fillSNP(char snp_str[], snp_t *snp_enc, int num_strains)
{

   
   /* fill in snp_str */

   
   for (int i = 0; i < num_strains; i++)
   {
      /* a filled zero is encoded as a 3, and a filled one is encoded as a 4:
           see fillin() */
      switch (snp_enc->iStatus)
      {
      case SNP_RIL_CONSTANT:
         if (snp_str[i] != toupper(snp_enc->cAlleleArray[0]))
         {
            snp_str[i] = tolower(snp_enc->cAlleleArray[0]);
         }
         break;
      /*
      case SNP_RIL_NOMATCH:
      case SNP_TOO_MANY_NT:
         snp_str[i] = '-';
         break;
      */
      case SNP_GOOD:
         if (snp_enc->iEncodedSNP[i] == 3)
         {
            snp_str[i] = tolower(snp_enc->cAlleleArray[0]);
            snp_enc->iEncodedSNP[i] = 2;
         }
         else if (snp_enc->iEncodedSNP[i] == 4)
         {
            snp_str[i] = tolower(snp_enc->cAlleleArray[1]);
            snp_enc->iEncodedSNP[i] = 2;
         }
         
         else if (snp_enc->iEncodedSNP[i] == 5)
         {
            snp_str[i] = 'n';
            snp_enc->iEncodedSNP[i] = 2;
         }
         break;
      default: 
         /* snp_str[i] = '-'; */
         break;
      }
   }
}


