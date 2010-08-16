/*
   File Name: prune.h
   System:
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
   Date developed:
   Purpose:
   Overview:
   Usage:
   Inputs:
   Returns:
   Effects:
   Assumptions:
   Dependencies:
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
 
 
#ifndef PRUNE_H
#define PRUNE_H

int runPruningIterations(em_params_t *em_params, snp_t *snps,
                         int emission_types, int prune_rule, int graph_flag, 
                         double prune_option, FILE *summary_file_fs);
 
int prune(int **prune_status_matrix, double **marg_prob_matrix, 
          snp_t *snps, int num_snps, int max_haplotypes, 
          double prune_rule);
           
void pruneLambda(int **prune_status_matrix, double ***transition_matrix, 
                 int num_snps, int max_haplotypes, int normalize);
void pruneOutp(int **prune_status_matrix, double ***emission_matrix, 
               int num_snps, int max_haplotypes, int emission_types);
               
void rowNormalize(double **p, int height, int width);

int **initPruneStatus(int num_snps, int max_haplotypes);
 
#endif
