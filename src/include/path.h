/*
   File Name: path.h
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
      Algorithm designed by Jin Szatkeiwicz under the supervision of 
      Gary Churchill; and a matlab version of the algorithm is available from 
      Jin.
   Date developed: July-October 2006
   Purpose: contains function prototypes for the path functions
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

#ifndef PATH_H
#define PATH_H

void maxSmthPath(const em_params_t *em_params, const em_out_t *em_out, 
                 const char *prefix, int smth_out, 
                 int **strain_paths, double **prob_matrix);
                 
void viterbiPath(const em_params_t *em_params, const em_out_t *em_out, 
                 int **strain_paths, double **log_prob_matrix, double **prob_matrix);
                 


#endif
