/*
   File Name: dataIO.h
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
      Algorithm designed by Jin Szatkeiwicz under the supervision of 
      Gary Churchill; and a matlab version of the algorithm is available from 
      Jin.
   Date developed: July-October 2006
   Purpose: function prototypes for externally referenced
      functions defined in dataIO.h
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


#ifndef DATAIO_H
#define DATAIO_H

typedef struct {

   int *iEncodedSNP;
   int  iNumHaplotypes;
   int  iStatus;
   int  index;
   char cAlleleArray[2];

} snp_t;

int     initializeDataSource(char *source, char *training_strains, 
                             double max_missingness, int num_parents, 
                             int warn_snp_removed, int filter_flags);

int count_good_snps(int chromosome, int filter_flags);

snp_t *loadSNPs(int chromosome, int max_haplotypes, int *good_snps, 
                int filter_flags);

double *calculateMissingness(snp_t *snps, int chromosome);

void    dumpDataInfo();

int     getNumStrains();
int     getChromosomeSize(int chromosome);
char  **getStrainNames();
int     getNumParents();
int     getNumTrainingStrains();

int     isTraining(int strain_index);



void writeLambdaOutp(char path[], double ***transition_matrix, 
                     double ***emission_matrix, double **marg_probability_matrix, 
                     int num_haplotypes, int num_snps, int miss);
                     
void writeFilledData(snp_t snps[], int chromosome, 
                     double **fill_probs, const char *prefix, 
                     int path_option, int sort_option, double threshold, 
                     int good_snps);
                     
void writePath(snp_t *snps, int **state_path, int num_snps, 
               int num_strains, const char prefix[], int path_option, 
               int chromosome, int sort_option, double ***emission_matrix);
               
void readLamboutpFile(char path[], double ***transition_matrix, 
                      double ***emission_matrix, int max_haplotypes, 
                      int emission_types, int num_snps);
                      
void writeLambdaOutpBinary(char path[], double ***transision_matrix, 
                           double ***emission_matrix,  
                           int max_haplotypes, int num_snps, int miss);
                           
void readLambdaOutpBinary(char path[], double ***transition_matrix, 
                          double ***emission_matrix, int max_haplotypes, 
                          int num_snps, int miss);
#endif
