/*
   File Name: state_sorting.c
   System: hmmSNP
   Programmer: Glen Beane
   Date developed:
   Purpose: provides state sorting functionality to hmmSNP
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
#include <sys/param.h>

#include "CLI.h"
#include "dataIO.h"
#include "sort_max_trace.h"
#include "hmm_em.h"
#include "constants.h"
#include "hmmUtil.h"
#include "graph.h"
#include "hmmMatrix.h"
#include "fill.h"

/*******************************************************************************
 * sort_states - wraps state sorting code, writing out *_sorted output files
 *
 * INPUTS:  snps - encoded snps
 *          num_good_snps - number of good snps
 *          em_params - parameters for hmm_em()
 *          em_output - output from hmm_em()
 *          transition_matrix - lambda
 *          emission_matrix 
 *          smth_out
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
void sort_states(snp_t *snps, int num_good_snps, em_params_t *em_params, 
                 em_out_t *em_output, double ***transition_matrix, 
                 double ***emission_matrix, int smth_out)
{ 
  char path[MAXPATHLEN];
  int num_haplotypes = get_num_haplotypes();
  int miss_option = get_miss_option();
  int num_emission_types;
  int path_option = get_path_option();
  double confidence_threshold = get_confidence_threshold();
  double **marginal_probability_matrix;
  double **fill_probability_matrix;
  int sort_option = get_sort_option();
  int chromosome = get_chromosome();
  
  /* right now SORT_MAX_TRACE is the only sorting algorithm implemented */
  assert(sort_option == SORT_MAX_TRACE);
  
  if (miss_option == MO_EMISSION) {
    /* miss_option is missing as emission type, so we have 3 possible emission
       types */
    num_emission_types = 3;
  }
  else {
    num_emission_types = 2;
  }

  /* sort based on the maximum trace */
  sort_max_trace(snps, em_output->last_transition_matrix, 
                 em_output->last_emission_matrix, num_good_snps, num_haplotypes, 
                 miss_option);


  /* compute marginal state probability */
  marginal_probability_matrix = allocate2Dd(num_good_snps, num_haplotypes);

  if (marginal_probability_matrix == NULL) {
    fprintf(stderr, 
    "    Unable to allocate memory for Marginal Probability Matrix\n");
    exit(EHMM_NOMEM);
  }  

  fill_probability_matrix = allocate2Dd(num_good_snps, em_params->num_strains);
  
  if (fill_probability_matrix == NULL) {
    fprintf(stderr, 
    "    Unable to allocate memory for Fill Probability Matrix\n");
    exit(EHMM_NOMEM);
  }

  computeMargProb(marginal_probability_matrix, em_output->last_transition_matrix, 
                  num_haplotypes, num_good_snps);
  
  snprintf(path, MAXPATHLEN-1, "%slamboutp_sorted.csv", em_params->file_prefix);
  
  writeLambdaOutp(path, em_output->last_transition_matrix, 
                  em_output->last_emission_matrix, marginal_probability_matrix, 
                  num_haplotypes, num_good_snps, num_emission_types);


  if (graphviz()) {
   snprintf(path, MAXPATHLEN-1, "%sgraph_sorted.dot", em_params->file_prefix);
   writeGraph(em_output->last_transition_matrix, em_output->last_emission_matrix, 
        marginal_probability_matrix, num_good_snps, num_haplotypes, 
        snps, num_emission_types, path, 100, 5);
  } 
  
  free2Dd(marginal_probability_matrix, num_good_snps);
        
              
  /* compute the simple path if necessary */
  if (path_option == PATH_MAXSMTH || path_option == PATH_BOTH) {
    fill_and_save(snps, em_params, em_output, PATH_MAXSMTH, 
                  sort_option, confidence_threshold,
                  num_emission_types, chromosome, smth_out,
                  fill_probability_matrix);
  }
   
   
  /* compute viterbi path if necessary (this is the most common path algorithm
     and is the default behavior */
  if (path_option == PATH_VITERBI || path_option == PATH_BOTH) {
    fill_and_save(snps, em_params, em_output, PATH_VITERBI,
                  sort_option, confidence_threshold,
                  num_emission_types, chromosome, smth_out,
                  fill_probability_matrix);
  }
  
  free2Dd(fill_probability_matrix, num_good_snps);
}
