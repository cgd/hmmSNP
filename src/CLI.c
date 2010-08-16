/*
   File Name: CLI.c
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
      Algorithm designed by Jin Szatkeiwicz under the supervision of 
      Gary Churchill; and a matlab version of the algorithm is available from 
      Jin.
   Date developed: July-October 2006
   Purpose: defines functions related to command line interface
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

#include <ctype.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <sys/param.h>

#include "constants.h"
#include "CLI.h"
#include "version.h"



static cliParams_t cliParams;


static int no_fix_outp_flag = FALSE;
static int prune1_flag    = FALSE;
static int prune2_flag    = FALSE;



/*******************************************************************************
 * cliParams_t *parseCL(int argc, char **argv)
 *  parse command line options and return parameters in cliParams_t struct
 *
 * INPUTS: int argc - argc passed to main function
 *         char **argv - argv passed to main
 *
 * RETURNS: pointer to cliParams_t struct containing parameters
 *
 * ASSUMES:
 *
 * EFFECTS: allocates memory for cliParams_t struct
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
cliParams_t *parseCL(int argc, char **argv)
{


   static struct option longopts[] = {
      {"annotation_columns", required_argument, NULL, 'C'},
      {"chromosome", required_argument, NULL, 'c'},
      {"emission_prior", required_argument, NULL, 'e'},
      {"graphviz", no_argument, NULL, 'g'},
      {"haplotypes", required_argument, NULL, 'h'},
      {"help", no_argument, NULL, 'u'},
      {"iterations", required_argument, NULL, 'i'},
      {"lambda_ratio", required_argument, NULL, 'l'},
      {"lambda_outp_file", required_argument, NULL, 'F'},
      {"keep_partial", no_argument, &(cliParams.keep_partial), TRUE},
      {"maxmiss", required_argument, NULL, 'm'},
      {"missopt", required_argument, NULL, 'o'},
      {"no_em", no_argument, &(cliParams.no_em), TRUE},
      {"parents", required_argument, NULL, 'p'},
      {"partial_output", required_argument, NULL, 'K'},
      {"path", required_argument, NULL, 'P'},
      {"prefix", required_argument, NULL, 'f'},
      {"prune1", optional_argument, NULL, 'y'},
      {"prune2", optional_argument, NULL, 'z'},
      {"pseudo_mod", required_argument, NULL, 'M'},
      {"pseudo_option", required_argument, NULL, 'O'},
      {"random_starts", required_argument, NULL, 'R'},
      {"random_iterations", required_argument, NULL, 'I'},
      {"sort", optional_argument, NULL, 'S'},
      {"sort_only", optional_argument, NULL, 'T'},
      {"stop", required_argument, NULL, 's'},
      {"confidence_threshold", required_argument, NULL, 'H'},
      {"tolerance", required_argument, NULL, 't'},
      {"tolerance_relaxed", required_argument, NULL, 'r'},
      {"training_strains", required_argument, NULL, 'Y'}, 
      {"version", no_argument, NULL, 'v'},
      {"remove_constant_snps", no_argument, &(cliParams.keep_constant), FALSE},
      {"no_fix_outp", no_argument, &no_fix_outp_flag, TRUE},
      {"smthout", no_argument, &(cliParams.smth_out), TRUE},
      {"warn_snp_removed", no_argument, &(cliParams.warn_if_removed), TRUE},
      {NULL, 0, NULL, 0}
   };

   int ch;
   int err = 0;
   char *end_ptr;
   char *token;
   char c;
   int i;
   int warn_hap_set = FALSE;
   int rand_starts_specified = FALSE;
   int emission_prior_read = 0;

   /* setup default values */
   cliParams.miss_option = MO_RANDOM;
   cliParams.max_iterations = MAX_ITERATIONS;
   cliParams.smth_out = 0;
   cliParams.num_haplotypes = 0;
   cliParams.warn_if_removed = FALSE;
   cliParams.keep_constant = TRUE;
   cliParams.fix_emission = FALSE;
   cliParams.random_starts = RANDOM_STARTS_DEFAULT;
   cliParams.random_iterations = RANDOM_ITERATIONS_DEFAULT;
   cliParams.sort_option = NO_SORT;
   cliParams.sort_only = FALSE;
   cliParams.graphviz = FALSE;
   cliParams.num_annotation_col = -1;
   cliParams.no_em = FALSE;
   cliParams.path_option = PATH_VITERBI;
   cliParams.stop_option = SO_CONV;
   cliParams.partial_output_interval = 0;
   cliParams.keep_partial = FALSE;
   cliParams.user_prefix_flag = FALSE;
   
   cliParams.emission_prior[0] = 9;
   cliParams.emission_prior[1] = 9;
   cliParams.emission_prior[2] = 1;

   cliParams.max_miss_rate = 0.8;
   cliParams.pseudo_option = PSEUDO_OPTION;
   cliParams.pseudo_mod = PSEUDO_MOD;
   cliParams.tolerance = TOLERANCE_DEFAULT;
   cliParams.tolerance_relaxed = TOLERANCE_RELAXED_DEFAULT;
   cliParams.lambda_ratio = LAMBDA_RATIO_DEFAULT;
   cliParams.confidence_threshold = CONFIDENCE_THRESHOLD_DEFAULT;
   
   cliParams.num_parent_strains = 0;
   
   cliParams.prune_rule = PRUNE_NONE;
   cliParams.prune_option = 0.0;

   
   cliParams.lambda_outp_file = NULL;
   cliParams.file_prefix = NULL;
   
   cliParams.training_strains = NULL;
   cliParams.training_strains_flag = FALSE;
   
   
   while ((ch = getopt_long(argc, argv, 
                            "c:e:gm:o:p:h:l:s:t:r:i:f:uvy:z:C:I:F:H:K:M:O:P:R:S:T:Y:", 
                            longopts, NULL)) != -1)
   {
   
      printf("\n");
      
      switch(ch)
      {
      case 'c':
      
         
         cliParams.chromosome = -1;
         if (toupper(optarg[0]) == 'X')
         {
            cliParams.chromosome = CHR_X;
         }
         else if (toupper(optarg[0]) == 'Y')
         {
            cliParams.chromosome = CHR_Y;
         }
         else if (toupper(optarg[0]) == 'M')
         {
            cliParams.chromosome = CHR_M;
         }
         else
         {
            cliParams.chromosome = (int)strtol(optarg, &end_ptr, 10);
            if (*end_ptr != '\0')
            {
               /* argument contained some invalid characters */
              printf(" * invalid chromosome\n");
              exit(EHMM_USAGE);
            }
         }
         
         if (cliParams.chromosome < 1 || cliParams.chromosome > CHR_MAX)
         {
            printf(" * invalid chromosome\n");
            err++;
         }

         break;
         
      case 'e':
         token = strtok(optarg, ",");
         i = 0;
         while (i < 3 && token != NULL)
         {
            cliParams.emission_prior[i] = (int)strtol(token, &end_ptr, 10);
            if (*end_ptr != '\0')
            {
               /*arg contained some invalid characters */
               printf(" * Invalid argument for emission_prior.\n");
               err++;
            }
            if (cliParams.emission_prior[i] < 1)
            {
               printf(" * Invalid argument for emission_prior: %d\n", 
                      cliParams.emission_prior[i]);
               printf("      each value must be >= 1\n");
               err++;
            }
            token = strtok(NULL, ",");
            i++;

         }
         emission_prior_read = i;
         break;
         
      case 'g':
         cliParams.graphviz = TRUE;
         break;
         
      case 'l':
         cliParams.lambda_ratio = strtod(optarg, &end_ptr);
         if (*end_ptr != '\0')
         {
            /* argument contained some invalid characters */
            printf(" * Invalid argument for lambda ratio.\n");
            err++;
         }
         else if (cliParams.lambda_ratio <= 0.0 
                  || cliParams.lambda_ratio >= 1.0)
         {
            printf(" * Lambda Ratio must be between zero and 1\n");
            err++;
         }
         break;
         
      case 'm':
         
         cliParams.max_miss_rate = strtod(optarg, &end_ptr);
         if (*end_ptr != '\0')
         {
            /* argument contained some invalid characters */
            printf(" * Invalid argument for Max miss rate.\n");
            err++;
         }  
         
         break; 

      case 'o':
         
         c = toupper(optarg[0]);
         if (c == 'C')
         {
            cliParams.miss_option = MO_COMPLETE;   
         }
         else if (c == 'R')
         {
            cliParams.miss_option = MO_RANDOM;                  
         }
         else if (c == 'E')
         {
            cliParams.miss_option = MO_EMISSION;                  
         }
         
         break;      

      case 'p':
         
         cliParams.num_parent_strains = (int)strtol(optarg, &end_ptr, 10);
         if (*end_ptr != '\0')
         {
            /* argument contained some invalid characters */
            printf(" * invalid argument for number of parent strains\n");
            err++;
         }
         if (no_fix_outp_flag == FALSE)
         {
            cliParams.fix_emission = TRUE;
         }
         if (cliParams.num_parent_strains != 2)
         {
            printf(" * Invalid number of parent strains: %d\n", 
               cliParams.num_parent_strains);
         }
         if (cliParams.num_haplotypes != 0) 
         {
            warn_hap_set = TRUE;
         }
         cliParams.num_haplotypes = cliParams.num_parent_strains;
         break;
         
      case 'h':
         if (cliParams.num_parent_strains >= 2)
         {
            warn_hap_set = TRUE;
         }
         else
         {
            cliParams.num_haplotypes = (int)strtol(optarg, &end_ptr, 10);
            if (*end_ptr != '\0')
            {
               /* argument contained some invalid characters */
               printf(" * invalid argument for haplotypes\n");
               err++;
            }
         }
         break;
         
      case 'r':
         printf("WARNING: pruning currently disabled, ignoring tolerance_relaxed option\n");
         break;
#if 0
         cliParams.tolerance_relaxed = strtod(optarg, &end_ptr);
         if (*end_ptr != '\0')
         {
            /* argument contained some invalid characters */
            printf(" * invalid argument for tolerance_relaxed\n");
            err++;
         }
         
         break;
#endif

      case 's':
         
         c = toupper(optarg[0]);
         if (c == 'I')
         {
            cliParams.stop_option = SO_MAXI;   
         }
         else if (c == 'C')
         {
            cliParams.stop_option = SO_CONV;                  
         }
         else
         {
            printf(" * invalid argument for --stop\n");
            err++;
         }
         
         if (optarg[1] != '\0')
         {
           /* should we print an error or warning here?  
              for now we'll just ignore anything after c or s*/
         }
         
         break;
         
      case 't':
         
         cliParams.tolerance = strtod(optarg, &end_ptr);
         if (*end_ptr != '\0')
         {
            /* argument contained some invalid characters */
            printf(" * invalid argument for tolerance\n");
            err++;
         }
         
         break;
         
      case 'i':
         
         cliParams.max_iterations = (int)strtol(optarg, &end_ptr, 10);
         if (*end_ptr != '\0')
         {
            /* argument contained some invalid characters */
            printf(" * Invalid argument for max iterations\n");
            err++;
         }
         
         break;
         
      case 'f':
         
         cliParams.file_prefix = optarg;
         cliParams.user_prefix_flag = TRUE;
         break;
         
      case 'u':
         
         usage();
         exit(0);
         
         break;
      case 'v':
         
         version();
         exit(0);
         
         break;
         
      case 'y':
         printf("WARNING:  state pruning is currently disabled, ignoring pruning option\n");
         break;
#if 0         
         if (cliParams.prune_rule != PRUNE_NONE)
         {
            /* warn previous prune rule value is getting over written */
         }
         cliParams.prune_rule = PRUNE_RULE_1;
         if (optarg != NULL)
         {
            cliParams.prune_option = strtod(optarg, &end_ptr);
            if (*end_ptr != '\0')
            {
               /* argument contained some invalid characters */
               printf(" * invalid argument for tolerance\n");
               err++;
            }
         }
         else
         {
            cliParams.prune_option = PRUNE_OPTION_1_DEFAULT;
         }
         break;
#endif

      case 'z':
         printf("WARNING:  state pruning is currently disabled, ignoring pruning option\n");
         break;
#if 0         
         if (cliParams.prune_rule != PRUNE_NONE)
         {
            /* warn previous prune rule value is getting over written */
         }
         cliParams.prune_rule = PRUNE_RULE_2;
         if (optarg != NULL)
         {
            cliParams.prune_option = strtod(optarg, &end_ptr);
            if (*end_ptr != '\0')
            {
               /* argument contained some invalid characters */
               printf(" * invalid argument for prune rule\n");
               err++;
            }
         }
         else
         {
            cliParams.prune_option = PRUNE_OPTION_2_DEFAULT;
         }      
         break;
#endif

      case 'C':
         cliParams.num_annotation_col = (int)strtol(optarg, &end_ptr, 10);
         if (*end_ptr != '\0' || cliParams.num_annotation_col < 1)
         {
            printf(" * Invalid argument for --annotation_columns.\n");
            err++;
         }
         break;
         
      case 'I':
         
         cliParams.random_iterations = (int)strtol(optarg, &end_ptr, 10);
         if (*end_ptr != '\0' || cliParams.random_iterations < 1)
         {
            /* argument contained some invalid characters */
            printf(" * Invalid argument for random_iterations\n");
            err++;
         }
         
         break;
         
      case 'F':
         cliParams.lambda_outp_file = optarg;
         break;
         
      case 'H':
         cliParams.confidence_threshold = strtod(optarg, &end_ptr);
         if (*end_ptr != '\0' || cliParams.confidence_threshold < 0
             || cliParams.confidence_threshold > 1)
         {
            /* argument contained some invalid characters */
            printf(" * Invalid argument for confidence_threshold\n");
            err++;
         }
      case 'K':
         cliParams.partial_output_interval = (int)strtol(optarg, &end_ptr, 10);
         if (*end_ptr != '\0' || cliParams.partial_output_interval < 0)
         {
            printf(" * Invalid argument for --partial_output.\n");
            err++;
         }
         break;
         
      case 'M':
         cliParams.pseudo_mod = strtod(optarg, &end_ptr);
          if (*end_ptr != '\0')
          {
             printf(" * Invalid argument for pseudo_mod\n");
             err++;
          }
         break;
         
      case 'O':
          cliParams.pseudo_option = strtod(optarg, &end_ptr);
          if (*end_ptr != '\0')
          {
             printf(" * Invalid argument for pseudo_option\n");
             err++;
          }
         break;
         
      case 'P':
         c = tolower(optarg[0]);
         if (c == 'b')
         {
            cliParams.path_option = PATH_BOTH;
         }
         else if (c == 's')
         {
            cliParams.path_option = PATH_MAXSMTH;
         }
         else if (c == 'v')
         {
            cliParams.path_option = PATH_VITERBI;
         }
         else
         {
            printf(" * invalid argument for --path\n");
            err++;
         }
         
         break;
         
      case 'R':
         cliParams.random_starts = (int)strtol(optarg, &end_ptr, 10);
         if (*end_ptr != '\0')
         {
            printf(" * Invalid argument for --random_starts.\n");
            err++;
         }
         rand_starts_specified = TRUE;
         break;
         
      case 'T':
         cliParams.sort_only = TRUE;
         /*fall through is intentional */
      case 'S':
         if (optarg != NULL)
         {
            /* currently the only supported option is max_trace */
            if (strcmp(optarg, "max_trace") == 0)
            {
               cliParams.sort_option = SORT_MAX_TRACE;
            }
            else
            {
               printf(" * Invalid sorting option.\n");
               err++;
            }
         }
         else
         {
            /* default */
            cliParams.sort_option = SORT_MAX_TRACE;
         }
         break;
         
      case 'Y':
         cliParams.training_strains = optarg;
         cliParams.training_strains_flag = TRUE;
         break;
         
      case 0:
         break;
      default:
         err++;
         break;
      }
      
   
   }
      
   if (emission_prior_read != 0) 
   {
      if (cliParams.miss_option == MO_RANDOM && emission_prior_read != 2)
      {
         printf("ERROR: \"Missing as random\" requires two emission_prior values, but %d were given.\n", 
                emission_prior_read);
         err++;
      }
      else if (cliParams.miss_option == MO_EMISSION && emission_prior_read != 3)
      {
         printf("ERROR: \"Missing as emission\" requires three emission_prior values, but %d were given.\n", 
                emission_prior_read);
         err++;
      }
   }
   
   if (prune1_flag)
   {
      cliParams.prune_rule = PRUNE_RULE_1;
      if (cliParams.prune_option == 0.0)
      {
         cliParams.prune_option = PRUNE_OPTION_1_DEFAULT;
      }
   }
   else if (prune2_flag)
   {
      cliParams.prune_rule = PRUNE_RULE_2;
      if (cliParams.prune_option == 0.0)
      {
         cliParams.prune_option = PRUNE_OPTION_2_DEFAULT;      
      }
   }
   
   cliParams.path = argv[optind];
   
   if (cliParams.path == NULL)
   {
      printf("ERROR: no dataset specified!\n");
      err++;
   }
   
   if (err > 0)
   {
      usage();
      exit(EHMM_USAGE);
   }
   
   if (warn_hap_set == TRUE)
   {
      printf("WARNING: number of haplotypes and number of parents set.\n" 
             "         ignoring -haplotypes and setting max number of "
             "haplotypes to number of parents.\n");
   }
   
   if (cliParams.num_parent_strains != 0)
   {
      cliParams.random_starts = 1;
      cliParams.keep_constant = FALSE;
   
      if (rand_starts_specified == TRUE)
      {
         printf("WARNING: RIL data, --random_starts ignored (random_starts "
                "fixed at 1)\n");
      }
   }
   
   if (cliParams.partial_output_interval == 0 && cliParams.keep_partial)
   {
      printf("WARNING: --keep_partial is set without --partial_output\n");
   }
   
   if (cliParams.num_haplotypes == 0)
   {
      cliParams.num_haplotypes = NUM_HAPLOTYPES_DEFAULT;
   }
   
   
   return &cliParams;

}

/*******************************************************************************
 * void dumpCLParams(void) - dump options to stdout
 *
 * INPUTS: void
 *
 * RETURNS: void
 *
 * ASSUMES:
 *
 * EFFECTS: prints parameter values to stdout
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void dumpCLParams(void)
{

   
   if (cliParams.chromosome == CHR_X)
      printf(" Chromosome X\n");
   else if (cliParams.chromosome == CHR_Y)
      printf(" Chromosome Y\n");
   else if (cliParams.chromosome == CHR_M)
      printf(" Chromosome M (mitochondria)\n");
   else
      printf(" Chromosome %d\n", cliParams.chromosome);
   
   printf(" MaxMissRate: %f\n", cliParams.max_miss_rate);
   
   printf(" MissOption: ");
   switch (cliParams.miss_option)
   {
      case MO_COMPLETE:
         printf("COMPLETE\n");
         break;
      case MO_RANDOM:
         printf("RANDOM\n");
         break;
      case MO_EMISSION:
         printf("EMISSION\n");
         
   }
   
   if (cliParams.miss_option == MO_RANDOM)
   {
      printf(" Emission Prior: %d,%d\n", cliParams.emission_prior[0], 
             cliParams.emission_prior[1]);
   }
   else if (cliParams.miss_option == MO_EMISSION)
   {
      printf(" Emission Prior: %d,%d,%d\n", cliParams.emission_prior[0], 
             cliParams.emission_prior[1], cliParams.emission_prior[2]);
   }
   
   printf(" Pseudo Option: %f\n", cliParams.pseudo_option);
   
   printf(" Pseudo Mod: %f\n", cliParams.pseudo_mod);
   
   printf(" Lambda Ratio: %f\n", cliParams.lambda_ratio);
   
   printf(" NumHaplotypes: %d\n", cliParams.num_haplotypes);
   
   printf(" StopOption: ");
   switch (cliParams.stop_option)
   {
      case SO_MAXI:
         printf("Maximum Iterations\n");
         break;
      case SO_CONV:
         printf("Convergence (or max iterations, whichever first)\n");
         break;
   }
   
   printf(" Tolerance: %f\n", cliParams.tolerance);
   
   printf(" MaxIterations: %d\n", cliParams.max_iterations);
   
   printf(" Random Starts (parameter initialization): %d\n", 
          cliParams.random_starts);
   printf(" Random Starts #iterations: %d\n", cliParams.random_iterations);
   
   printf(" PathOption: ");
   switch (cliParams.path_option)
   {
      case PATH_MAXSMTH:
         printf("Simple path\n");
         break;
      case PATH_VITERBI:
         printf("Viterbi path\n");
         break;
      case PATH_BOTH:
         printf("Viterbi and Simple paths\n");
         break;
   
   }
   
   printf(" SmthOutOption: %d\n", cliParams.smth_out);
   
   printf(" Path: %s\n", cliParams.path);
   
   if (cliParams.file_prefix != NULL)
   {
      printf(" OutputFile Prefix: %s\n", cliParams.file_prefix);
   }
   
/*
   if (cliParams.prune_rule != PRUNE_NONE)
   {
      printf(" Prune Rule:  rule %d\n", cliParams.prune_rule);
   }
   else
   {
      printf(" Prune Rule: NONE\n");
   }
*/
   printf(" Fix Outp: ");
   (cliParams.fix_emission) ? printf("TRUE\n") : printf("FALSE\n");
   
   
   printf(" Sorting Option: ");
   switch (cliParams.sort_option)
   {
      case SORT_MAX_TRACE:
         printf("max_trace\n");
         break;
      case NO_SORT:
         printf("no sorting\n");
         break;
   }
   
   printf(" Keep constant SNPs: ");
   cliParams.keep_constant ? printf("TRUE\n") : printf("FALSE\n");\
   
   if (cliParams.training_strains != NULL)
   {
      printf(" Training columns: %s\n", cliParams.training_strains);
   }
 
   printf(" Graphviz output: ");
   (cliParams.graphviz) ? printf("TRUE\n") : printf("FALSE\n");
   
   printf(" Confidence Threshold: %f\n", cliParams.confidence_threshold);
   
   printf(" Partial Output Interval: %d\n", cliParams.partial_output_interval);
   
   printf(" Keep Partial Output: ");
   cliParams.keep_partial ? printf("TRUE\n") : printf("FALSE\n");
   
   printf("\n");
}

/*******************************************************************************
 * void usage(void)
 *
 * INPUTS: void
 *
 * RETURNS: void
 *
 * ASSUMES: none
 *
 * EFFECTS: prints helpful message to stdout
 *
 * ERROR CONDITIONS: none
 *
 * COMMENTS: this uses one printf statement and passes a concatenation of a 
 *    bunch of different strings.  hopefully we don't run into any type of 
 *    system limit of the size of a string that printf will accept as input
 *
 *    should remove hard coded default values and use values from constants.h
 *
 ******************************************************************************/
void usage()
{
   printf(
      "\n\nhmmSNP [options] [input_data]\n"
      "   Options:\n"
      "     --version (-v):       print version information and exit\n"
      "     --help (-u):          help / useage info.  Prints this message.\n\n"    
      "     --chromosome (-c):    specify chromosome to process (1-19,x,y,m). "
      "Required\n"
      "     --annotation_columns: specify number of annotation columns for "
      "single .csv\n"
      "                            input file\n"
      "     --haplotypes (-h):    number of haplotypes. (default 3)\n"
      "     --parents (-p):       number of parent strains (for RIL data)\n"
      "                             currently only --parents 2 is supported\n"
      "                             --parents 2 implies --haplotypes 2\n"
      "     --no_fix_outp:        if --parents has been passed, outp will "
      "be fixed.\n"
      "                             this option allows outp to be updated during\n"
      "                             EM like non RIL data sets\n"
      "     --training_strains    comma delimited list of strains to use for\n"
      "                             training, remaining strains will be used\n"
      "                             as test data. ranges are allowed (e.g. 1-5)\n"
      "                             strain columns are numbered from 1 to #strains\n");
   printf(
      "     --maxmiss (-m):       maximum allowable percentage of missing\n"
      "                             nucleotides, default=0.8 (80%%)\n"
      "     --path (-P):          path algorithm. (v)iterbi,(s)imple, or (b)oth\n"
      "                             (default v)\n"
      "     --smthout:            print final smth values (simple path only)\n\n"
      "     --prefix (-f):        override default output file prefix\n"
      "     --confidence_threshold:\n"
      "                           remove imputed genotypes with a confidence score\n"
      "                             less than this value from the filled_filtered.csv\n"
      "                             file. (default 0.6)\n"
      "     --graphviz (-g):      create graphviz output file\n\n"
      "     --missopt (-o):       missing option\n" 
      "                             -o c (complete data)\n"
      "                             -o r (random missing, default)\n"
      "                             -o e (missing as emission)\n"
      "     --emission_prior (-e):\n"
      "                           emission prior values. A comma delimited list\n"
      "                             of integer values (2 for missing as random,\n"
      "                             3 for missing as emission). Default: 9,9 for\n"
      "                             missing as random, 9,9,1 for missing as emission.\n"      
      "     --stop (-s):          stop option\n"
      "                             -s c (convergence OR max iterations, default)\n"
      "                             -s i (max iterations)\n");
   printf(
      "     --tolerance (-t):     convergence tolerance (default 0.000001 change in\n"
      "                             log-likelihood)\n"
#if 0
/* removed because pruning is currently broken */
      "     --tolerance_relaxed:  relaxed tolerance used for trimming models\n"
      "                             (default 0.0001)\n"
#endif
      "     --iterations (-i):    max iterations (default 10000)\n\n"         
      "     --random_starts:      specify number of random initializations\n"
      "                             (default 10, 1 for RIL data. can not be\n"
      "                             changed for RIL data)\n"
      "     --random_iterations:  specifies number of interations to use for\n"
      "                             random initializations (default 100).\n"
      "                             Ignored if random initialization is not done.\n\n"
      "     --lambda_ratio (-l):  specify lambda.ratio between 0 and 1 (default 0.96)\n"
      "                             for state transitions, the prior density of \n"
      "                             transition probabilities is biased towards the\n"
      "                             transitions between the same haplotype, with\n"
      "                             probability lambda, and is equally distributed \n"
      "                             among the other (H-1) haplotypes\n"
      "     --lambda_outp_file:   specify a csv file to use for initial values of lambda\n"
      "                             and outp (emission matrix) \n\n" 
      "     --pseudo_option: (-O) pseudocount. Default 0.1, which means adding 1/10th\n"
      "                             of sequence.\n");
   printf(
      "     --pseudo_mod (-M):    when contributing the pseudocount to the emission\n"
      "                             matrix we first multiply it by the pseudomod.\n"
      "                             (default 1.0)\n\n"
#if 0 
/* pruning removed, currently broken */
      "     --prune1:             trim using prune rule 1\n"
      "            [=optional rule]\n"
      "     --prune2:             trim using prune rule 2 (option/#strains)\n"
      "            [=optional rule]\n\n"
#endif
      "     --warn_snp_removed:   print warning message if a SNP is removed\n"
      "     --remove_constant_snps: remove any SNP whose ovserved genotypes are same\n" 
      "                               among all sequences. SNPs are removed before\n"
      "                               running HMM (note: constant SNPs are always\n"
      "                               removed with RIL data)\n\n"
      "     --sort:               sort states using default sorting" 
      "(max trace)\n"
      "            [=optional sorting algorith]\n"
      "     --sort_only:          do not run HMM, read in lambda/outp and sort.\n"
      "                [=optional sorting algorithm]\n"
      "                            default algorithm = max_trace\n"
      "                            requires --lambda_outp_file\n\n"
      "     --no_em:             compute path and fill in using parameters passed\n"
      "                            with --lambda_outp_file (do not run EM training\n\n"
      "     --partial_output:    specify number of interations between writing out\n"
      "                            partial results (default 0, no partial results).\n"
      "                            Previous partial results files will be deleted after\n"
      "                            a new file is sucessfully written to disk or\n"
      "                            after final results are written to disk.\n"
      "     --keep_partial:      keep partial output files (see above)\n");
    
}




/* the following functions are "getters" that were added so that we didn't have 
   to change function  parameters all over the place when adding a couple 
   new options in version 1.1.0, I think it would be good to add more of these
   rather than passing around a pointer to the cliParams struct */

/*******************************************************************************
 * int get_snp_offset()
 *
 * INPUTS:  void
 *
 * RETURNS: offset of SNP column
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
int get_snp_offset(void)
{
   if (cliParams.num_annotation_col != -1)
   {
      return cliParams.num_annotation_col + 1;
   }

   return SNP_OFFSET;

}

/*******************************************************************************
 * int get_num_haplotypes()
 *
 * INPUTS:  void
 *
 * RETURNS: number of haplotypes specified
 *
 * ASSUMES: This function assumes that that cliParams has been initialized.
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
int get_num_haplotypes(void)
{
   return cliParams.num_haplotypes;
}

/*******************************************************************************
 * int get_miss_option()
 *
 * INPUTS: void
 *
 * RETURNS: integer number representing the miss option
 *
 * ASSUMES: this function assumes that cliParams has been initialized.
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
int get_miss_option(void)
{
   return cliParams.miss_option;
}

/*******************************************************************************
 * int get_path_option()
 *
 * INPUTS: void
 *
 * RETURNS: path option
 *
 * ASSUMES: cliParams has been initialized
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
int get_path_option(void)
{
   return cliParams.path_option;
}

/*******************************************************************************
 * int get_sort_option()
 *
 * INPUTS: void
 *
 * RETURNS: sort option (NO_SORT, SORT_MAX_TRACE)
 *
 * ASSUMES: cliParams has been initialized
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
int get_sort_option(void)
{
   return cliParams.sort_option;
}

/*******************************************************************************
 * int get_chromosome()
 *
 * INPUTS: void
 *
 * RETURNS: chromosome value (integer code)
 *
 * ASSUMES: cliParams has been initialized
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
int get_chromosome(void)
{
   return cliParams.chromosome;
}

/*******************************************************************************
 * int get_num_emission_types()
 *
 * INPUTS: void
 *
 * RETURNS: integer number of emission types (2 or 3)
 *
 * ASSUMES: that cliParams has been initialized
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
int get_num_emission_types(void)
{
   if (cliParams.miss_option == MO_EMISSION)
   {
      return 3;
   }
   
   return 2;
}

/*******************************************************************************
 *  double get_pseudo_mod(void)
 *
 * INPUTS: void
 *
 * RETURNS: pseudomod value
 *
 * ASSUMES: parseCL has been called once
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS: 
 *
 ******************************************************************************/
double get_pseudo_mod(void)
{
   return cliParams.pseudo_mod;
}


/*******************************************************************************
 *  void get_emission_prior(int *emission_prior)
 *
 * INPUTS: int *emission_prior - pointer to integer array to store emission prior
 *
 * RETURNS: void,  emission prior values returned by reference (emission_prior)
 *
 * ASSUMES: parseCL has been called once and that emission_prior can hold 3 ints
 *
 * EFFECTS: modifies contents of emission_prior
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS: 
 *
 ******************************************************************************/
void get_emission_prior(int *emission_prior)
{
   int i;
   for (i = 0; i < 3; i++)
   {
      emission_prior[i] = cliParams.emission_prior[i];
   }
}


/*******************************************************************************
 * double get_confidence_threshold()
 *
 * INPUTS: void
 *
 * RETURNS: double confidence_threshold
 *
 * ASSUMES: that cliParams has been initialized
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
double get_confidence_threshold(void)
{
   return cliParams.confidence_threshold;
}

/*******************************************************************************
 * int graphviz()
 *
 * INPUTS: void
 *
 * RETURNS: value of graphviz attribute (true if we are to produce graphviz .dot 
 *   for output.
 *
 * ASSUMES: cliParams has been initialized.
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
int graphviz(void)
{
  return cliParams.graphviz;
}