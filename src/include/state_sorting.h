/*
   File Name: state_sorting.h
   System: hmmSNP
   Programmer: Glen Beane
   Date developed: August 2009
   Purpose: header file for state_sorting
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

#ifndef STATE_SORTING_H
#define STATE_SORTING_H

void sort_states(snp_t *snps, int num_good_snps, em_params_t *em_params, 
                 em_out_t *em_out, double ***transition_matrix, 
                 double ***emission_matrix, int smth_out);

#endif