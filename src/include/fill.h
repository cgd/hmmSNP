/*
   File Name: fill.h
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
      Algorithm designed by Jin Szatkeiwicz under the supervision of 
      Gary Churchill; and a matlab version of the algorithm is available from 
      Jin.
   Date developed: July-October 2006
   Purpose: this file contain function prototypes for externally referenced 
     functions defined in fill.c
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

#ifndef FILL_H
#define FILL_H

#include "dataIO.h"
#include "hmm_em.h"


void fill_and_save(snp_t *snps, em_params_t *em_params, em_out_t *em_out, 
                   int path_option, int sort_option, double threshold,
                   int num_emission_types, int chromosome, int smooth_out,
                   double **fill_probability_matrix);
                   
void  fillin(snp_t *snps, int **state_paths, double ***emission_matrix, 
            int num_snps, int num_strains, int miss, double **fill_prob_matrix);
            
void calculate_fill_rule(snp_t snp, double **emission_for_snp, 
                         int haplotypes, int emission_types, char fill_rule[]);
            
void fillSNP(char snp_str[], snp_t *snp_ptr, int num_strains);

#endif
