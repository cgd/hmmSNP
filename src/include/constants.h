/* test 
   File Name: constants.h
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
      Algorithm designed by Jin Szatkeiwicz under the supervision of 
      Gary Churchill; and a matlab version of the algorithm is available from 
      Jin.
   Date developed: July-October 2006
   Purpose: defines constants used in hmm program.
   Overview: This file contains a number of C preprocessor (cpp) defines. 
      These are intended to be included in every other .c file.
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

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "config.h"



#define HMM_VERSION VERSION
#define HMM_RELEASE_DATE "2/3/2010"

/* parameter defaults */
#define NUM_HAPLOTYPES_DEFAULT 3

/* error codes */
#define EHMM             1000
#define EHMM_NODATA      1001
#define EHMM_MISSINGDATA 1002
#define EHMM_NOMEM       1003
#define EHMM_FORMAT      1004
#define EHMM_PARSE       1005
#define EHMM_USAGE       1006
#define EHMM_NOPARENT    1007
#define EHMM_NOGOODSNPS  1008
#define EHMM_NOPARAM     1009
#define EHMM_SORT        1010


/* stop options */
#define SO_MAXI 0
#define SO_CONV 1

/* path options */
#define PATH_MAXSMTH 0
#define PATH_VITERBI 1
#define PATH_BOTH    2

/* missing options */
#define MO_COMPLETE 0
#define MO_RANDOM   1
#define MO_EMISSION 2

/* pruning options */
#define PRUNE_NONE    0
#define PRUNE_RULE_1  1
#define PRUNE_RULE_2  2

#define PRUNE_OPTION_1_DEFAULT 0.05
#define PRUNE_OPTION_2_DEFAULT 0.5

/* sorting options */
#define NO_SORT        0
#define SORT_MAX_TRACE 1

/* this program is designed for mice, we expect a maximum of 22 different 
   values for the chromosome option: chr 1-19, X, Y, M (mitochondria) */
#define NUM_CHROMOSOMES 22

/* buffer sizes */
#define LINE_BUFFER_SIZE 64000
#define MAX_PREFIX 128
#define LIST_SIZE 256

/* parent column offsets, to aid in code readability */
#define PARENT_1 0
#define PARENT_2 1

/* default column of first sequence. */
#define SNP_OFFSET 5

#define MAX_ITERATIONS 10000

/* constants used to initialize transition and emission matrixes */
#define RANDOM_STARTS_DEFAULT 10
#define RANDOM_ITERATIONS_DEFAULT 100
#define LAMBDA_RATIO_DEFAULT 0.96

/* default tolerance for hmm_em() function */
#define TOLERANCE_DEFAULT 0.000001

/* default tolerance for "pruning" models */
#define TOLERANCE_RELAXED_DEFAULT 0.0001

/* default confidence threshold, when generating the filled_filtered output file
   any imputed genotypes with a confidence score below this threshold will be
   removed */
#define CONFIDENCE_THRESHOLD_DEFAULT 0.6

/* chromosome numbers for X, Y, M. to be used as indexes into a zero 
   based array (e.g. our array of chromosome sizes) you must subtract 1 */
#define CHR_X NUM_CHROMOSOMES - 2
#define CHR_Y NUM_CHROMOSOMES - 1
#define CHR_M NUM_CHROMOSOMES 
#define CHR_MAX NUM_CHROMOSOMES

/* use pseudoption * PSEUDO_MOD when computing ssy in EM algorithm */
#define PSEUDO_MOD 1
#define PSEUDO_OPTION 0.1

/* any edge in the graph with a weight < EDGE_THRESHOLD will be hidden */ 
#define EDGE_THRESHOLD 0.05

#define TRUE 1
#define FALSE 0

/* format strings */
#define SUMMARY_FORMAT "%6d   %9.2e   %11d  %11.2f   %16d   %16d   %15.4f   %15.4f   %15.4f   %15.4f\n"

/* print extra diagnostic stuff
   right now there are two levels, zero (no extra output) and anything else */
#define HMM_VERBOSITY 0

/* SNP status values */
#define SNP_GOOD                  0
#define SNP_TOO_MANY_MISSING      1
#define SNP_BAD_FORMAT            2
#define SNP_TOO_MANY_NT           3
#define SNP_ONE_NT                4
#define SNP_RIL_CONSTANT          5
#define SNP_RIL_MISSING_PARENT    6
#define SNP_RIL_MISSING_PARENTS   7
#define SNP_RIL_NOMATCH           8

#define KEEP_CONSTANT_FLAG        0x00000001


#define TMP_FILE_NAME_BUFFER      128

/* TODO, hmmSNP should respect $TMP_DIR if set */
#define DEFAULT_TEMP_DIR "/tmp"

#endif
