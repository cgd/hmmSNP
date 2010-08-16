/*
   File Name:
   System:
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
   Date developed:
   Purpose:
   Overview: This file contains various utility functions that don't necessarily
             belong in any of the other function groups
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

#ifndef HMMUTIL_H
#define HMMUTIL_H

#define EUSAGE 1001
#define ELSTAT 1002
 
/* string related functions */ 
int isDir(const char *path);
int getChrNum(char *chr_str);
void getSNP(const char *line, char snp_str[], int offset, 
            int strain_flags[]);

char *setupOutputPrefix(int chromosome, char prefix[]);


/* memory allocation and free functions */

/* ... for doubles */
double **allocate2Dd(int height, int width);
void free2Dd(double **p, int rows);
/* ... for integers */
int **allocate2Di(int height, int width);
void free2Di(int **p, int rows);


int fprint_chr(FILE *stream, int chr);

int expandStrainList(char *list, int num_strains_total, 
                     char **expanded_list, int *num_strains_expanded);

#endif
