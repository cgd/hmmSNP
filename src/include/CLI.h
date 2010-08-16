/*
   File Name: CLI.h
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
      Algorithm designed by Jin Szatkeiwicz under the supervision of 
      Gary Churchill; and a matlab version of the algorithm is available from 
      Jin.
   Date developed: July-October 2006
   Purpose: contains function related to command line interface
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


#ifndef CLI_H
#define CLI_H



typedef struct 
{

   int chromosome;
   int miss_option;
   int num_haplotypes;
   int stop_option;
   int max_iterations;
   int path_option;
   int smth_out;
   int num_parent_strains;
   int fix_emission;
   int warn_if_removed;
   int random_starts;
   int random_iterations;
   int keep_constant;
   int sort_option;
   int sort_only;
   int training_strains_flag;
   int graphviz;
   int num_annotation_col;
   int no_em;
   int partial_output_interval;
   int keep_partial;
   int user_prefix_flag;
   int emission_prior[3]; 
  
   double tolerance;
   double tolerance_relaxed;
   
   double max_miss_rate;
   double lambda_ratio;
   double pseudo_option;
   
   double pseudo_mod;
   
   double confidence_threshold;
   
   int    prune_rule;
   double prune_option;
   
   char *path;
   char *file_prefix;
   char *lambda_outp_file;
   char *training_strains;
   

} cliParams_t;


cliParams_t *parseCL(int argc, char **argv);

void dumpCLParams(void);

void usage(void);

void version(void);

int get_snp_offset(void);
int get_num_haplotypes(void);
int get_miss_option(void);
int get_path_option(void);
int get_sort_option(void);
int get_chromosome(void);
int get_num_emission_types(void);
double get_pseudo_mod(void);
double get_confidence_threshold(void);

void get_emission_prior(int *);
int graphviz();

#endif
