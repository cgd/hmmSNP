/*
   File Name: main.c
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
      Algorithm designed by Jin Szatkeiwicz under the supervision of 
      Gary Churchill; and a matlab version of the algorithm is available from 
      Jin.
   Date developed: July-October 2006
   Purpose: includes main function - major program structure and flow control
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
 
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>

#include "CLI.h"
#include "constants.h"
#include "dataFilter.h"
#include "dataIO.h"
#include "fill.h"
#include "filter.h"
#include "graph.h"
#include "hmm_em.h"
#include "hmmMatrix.h"
#include "hmmUtil.h"
#include "path.h"
#include "prune.h"
#include "state_sorting.h"
#include "temp_file_cleanup.h"




int main(int argc,  char **argv)
{

   /* struct to hold command line parameters */
   cliParams_t *cliParams;
   
   /* input and output from hmm_em */
   em_params_t emParams;
   em_out_t   *emOutput;
   
   double **marg_prob_matrix;
   
   /* state pruning option, default to no pruning */
   double prune_option = PRUNE_NONE;
   
   /* 2D array of confidence scores for filled alleles */
   double **fill_confidences;
   
   /* final transition and emission matrixes from hmm_em function */
   double ***last_transition_matrix;
   double ***last_emission_matrix;
   
   /* total number of SNPs for the specified chromosome in the input file */
   int num_SNPs;
   
   /* number of SNPs that pass our cleaning process */
   int num_good_snps;
   
   /* number of strains */
   int num_strains;
   
   /* number of emission types; either 2 (default), or 3 (missing as emission option) */
   int num_emission_types;
   
   int model_number = 0;

   /* number of pruning iterations to run, specified by user */
   int pruning_iterations;
   
   /* SNP filtering options */
   int filter_flags = 0;
   
   char path[MAXPATHLEN];
   char *file_prefix;
   
   /* SNP data - dynamicaly allocated array of our snp_t structs */
   snp_t *SNPs;


   FILE *fsSummaryFile;


   

   /* parse command-line options, store options in cliParams struct */
   cliParams = parseCL(argc, argv);
   
   /* if we are keeping constant SNPs, then set the appropriate filter flag */
   if (cliParams->keep_constant == TRUE)
   {
      filter_flags = filter_flags | KEEP_CONSTANT_FLAG;
   }
   
   /* initialize the datasource */
   initializeDataSource(cliParams->path, cliParams->training_strains, 
                        cliParams->max_miss_rate, 
                        cliParams->num_parent_strains, 
                        cliParams->warn_if_removed, filter_flags);
                        
   /* calculate the number of SNPs in the input file for this chromosome */
   num_SNPs = getChromosomeSize(cliParams->chromosome);
   
   
   num_strains = getNumStrains();
   
   /* print overview of command line options */
   printf("\n--------------------------------------------------------------\n");   
   dumpCLParams();
   printf("--------------------------------------------------------------\n"); 
   dumpDataInfo();
   printf("--------------------------------------------------------------\n");


   /* exit if there were no SNPs for the specified chromosome */
   if (num_SNPs == 0)
   {
      printf("No SNPs for specified chromosome found.\n");
      exit (EHMM_NOGOODSNPS);
   }


   /* calculate the prune_option if a prune rule was specified */
   if (cliParams->prune_rule == PRUNE_RULE_1)
   {
      prune_option = cliParams->prune_option;
   }
   else if (cliParams->prune_rule == PRUNE_RULE_2)
   {
      prune_option = cliParams->prune_option / num_strains;
   }
   
   
   /* read in data for chromosome of interest */
   SNPs = loadSNPs(cliParams->chromosome, cliParams->num_haplotypes, 
                   &num_good_snps, filter_flags);
   
   printf("  ** %d SNPs removed\n", num_SNPs - num_good_snps);
   
   if (num_good_snps == 0)
   {
      printf("No good SNPs found.\n");
      exit (EHMM_NOGOODSNPS);
   }

     
   /* generate a prefix for output files based on the chromosome number or 
      a user specified prefix */
   file_prefix = setupOutputPrefix(cliParams->chromosome, 
                                   cliParams->file_prefix);



   /* if miss option is COMPLETE make sure the data really is complete. */  
   if (cliParams->miss_option == MO_COMPLETE)
   {
      for (int i = 0; i < num_good_snps; i++)
      {
         if (getMissingness(&SNPs[i]) != 0.0)
         {
            printf("ERROR: Missing values found in data,\n but Missing Option " 
                   "\"Complete\" specified. (-missopt c).\n");
            exit(EHMM_MISSINGDATA);
         }
      }
   }
   

   /* open and setup summary file */
   snprintf(path, MAXPATHLEN-1, "%ssummary.txt", file_prefix);
   fsSummaryFile = fopen(path, "w");
   
   if (fsSummaryFile == NULL)
   {
      perror("Unable to open summary output file");
      exit(errno);
   }
   
   fprintf(fsSummaryFile, "hmmSNP VERSION %s SUMMARY FILE\n", HMM_VERSION);
   fprintf(fsSummaryFile, "Data File: %s\n", cliParams->path);
   fprintf(fsSummaryFile, "Chromosome: "); 
   fprint_chr(fsSummaryFile, cliParams->chromosome);
   fprintf(fsSummaryFile, "\nNumber of SNPs in chromosome: %d\n", num_SNPs);
   fprintf(fsSummaryFile, "Humber of Haplotypes: %d\n", 
           cliParams->num_haplotypes);
   fprintf(fsSummaryFile, "Model#   Tolerance   #iterations      time(s)"
                          "            #States            #Params        "
                          "likelihood               AIC              BIC1"
                          "              BIC2\n");
                          
   fflush(fsSummaryFile);

   
   /* determine number of emission types */
   num_emission_types = get_num_emission_types();   



   /* setup E-M input parameters */

   /* data parameters */
   emParams.snps = SNPs;
   emParams.num_snps = num_good_snps;
   emParams.num_strains = num_strains;


   /* parameters from user */
   emParams.max_haplotypes = cliParams->num_haplotypes;
   emParams.fix_emission = cliParams->fix_emission;
   emParams.miss_option = cliParams->miss_option;
   emParams.stop_option = cliParams->stop_option;
   emParams.partial_output_interval = 0;
   emParams.keep_partial = cliParams->keep_partial;
   
   /* hmm_em will allocate these if they are NULL */
   emParams.input_transition_matrix = NULL;
   emParams.input_emission_matrix = NULL;
   
   /* initialize all states to "unpruned" */
   emParams.prune_status_matrix = initPruneStatus(num_good_snps, 
                                                  cliParams->num_haplotypes);
   emParams.fix_emission = cliParams->fix_emission;
  
   emParams.lambda_ratio = cliParams->lambda_ratio;
   emParams.pseudo_option = cliParams->pseudo_option;
   emParams.file_prefix = file_prefix;
   
   /* do we need to run hmm_em?, if we aren't just sorting or filling with 
      existing transition and emission matrixes, then run hmm_em */
   if (cliParams->sort_only == FALSE && cliParams->no_em == FALSE)
   {
      /* setup emParams with the proper values */
   
   
      /* did the user provide initial values ? */
      if (cliParams->lambda_outp_file != NULL)
      {     
         /* need to allocate memory for lambda and outp */
         emParams.input_transition_matrix = allocateLambda(num_good_snps,  
                                                      cliParams->num_haplotypes);
         emParams.input_emission_matrix = allocateOutp(num_good_snps, 
                                                      cliParams->num_haplotypes,
                                                      num_emission_types);
         readLamboutpFile(cliParams->lambda_outp_file, 
                          emParams.input_transition_matrix, 
                          emParams.input_emission_matrix, 
                          cliParams->num_haplotypes, num_emission_types, 
                          num_good_snps);
      }
      else
      {
         printf("Setting up Initial Lambda and Output Matrixes...");
         emParams.max_iterations = cliParams->random_iterations;
         
         /* this next step can take a while. make sure the stdout buffer is 
            flushed, so the user can tell what they are waiting on */
         fflush(stdout);   
         randStart(&emParams, cliParams->random_starts, file_prefix); 
         printf("Done.\n");
      }   
   
      emParams.max_iterations = cliParams->max_iterations;
   
   
   
      /* do E-M */
      
      printf("Doing EM training...\n");

      /* pruning command line options are currently disabled, so this block 
         can't be executed right now */
      if (cliParams->prune_rule != PRUNE_NONE)
      {
         
         printf(" Starting Pruning Models...\n");
         
         emParams.tolerance = cliParams->tolerance_relaxed;
          
         pruning_iterations = runPruningIterations(&emParams,
                                                   SNPs, 
                                                   num_emission_types, 
                                                   cliParams->prune_rule,
                                                   cliParams->graphviz, 
                                                   cliParams->prune_option,
                                                   fsSummaryFile);
          
         /* 
           emParams.input_transition_matrix == emOut->last_transition_matrix
           emParams.input_emission_matrix == emOut->last_emission_matrix
           
           as of last call to hmm_em in runPruningIterations, so now the 
           Lambda and Outp input matrixes have been "pruned"
           
          */
      
         printf(" Starting final model...\n");
      
      }

      snprintf(path, MAXPATHLEN-1, "%sloglike.txt", file_prefix);
      emParams.tolerance = cliParams->tolerance;
      emParams.partial_output_interval = cliParams->partial_output_interval;

      /****** this is where the magic happens - run EM ******/
      emOutput = hmm_em(&emParams, path);


      
      /* check for convergence */
      if (emOutput->err == ELIKE)
      {
         printf("    did not converge: loglikelihood decreased at iteration %d\n", 
                emOutput->last_iteration);
      }
      else
      {
         printf("    converged at %d\n", emOutput->last_iteration);
      }
      
      /* print timing information */
      printf("    E-M Elapsed Time (in iterations): %.2f(s)\n", 
             emOutput->elapsed_time);
      printf("    Avg Iteration: %.3f(s)\n", 
             emOutput->elapsed_time / emOutput->last_iteration);
       

      fprintf(fsSummaryFile, SUMMARY_FORMAT, model_number, emParams.tolerance, 
              emOutput->last_iteration, emOutput->elapsed_time, 
              emOutput->num_states, emOutput->k, emOutput->last_log_liklihood, 
              emOutput->aic, emOutput->bic[0], emOutput->bic[1]);
      fflush(fsSummaryFile);

      marg_prob_matrix = allocate2Dd(num_good_snps, cliParams->num_haplotypes);

      if (marg_prob_matrix == NULL)
      {
         fprintf(stderr, 
              "    Unable to allocate memory for Marginal Probability Matrix\n");
         exit(EHMM_NOMEM);
      }  

      computeMargProb(marg_prob_matrix, emOutput->last_transition_matrix, 
                      cliParams->num_haplotypes, num_good_snps);
                     
      /* write results out to a file */
      
      /* setup filename */
      snprintf(path, MAXPATHLEN-1, "%slamboutp.csv", file_prefix);          

      /* write out the lambda & outp info */
      writeLambdaOutp(path, emOutput->last_transition_matrix, 
                      emOutput->last_emission_matrix, marg_prob_matrix, 
                      cliParams->num_haplotypes, emParams.num_snps, 
                      emOutput->miss);
                
      /* write graphviz file */   
      if (cliParams->graphviz == TRUE)
      {
         /* setup filename, if pruning include model number */
         printf("   Writing Graphviz file...");
         snprintf(path, MAXPATHLEN-1, "%sgraph.dot", file_prefix);
         
         
         writeGraph(emOutput->last_transition_matrix, emOutput->last_emission_matrix, 
                    marg_prob_matrix, num_good_snps, cliParams->num_haplotypes, 
                    SNPs, emOutput->miss, path, 100, 5);
         printf("Done.\n");
      }

   
      free2Dd(marg_prob_matrix, num_good_snps);
      last_transition_matrix = emOutput->last_transition_matrix;
      last_emission_matrix = emOutput->last_emission_matrix;
   }
   else /*  cliParams->sort_only == TRUE  or cliParams->no_em == TRUE*/
   {
      double  **dPredMatrix;
      double  **dFiltMatrix;
      double  **dSmthMatrix;
      
      if (cliParams->lambda_outp_file != NULL)
      {
         /* read in lambda/outp values */
      
         /* need to allocate memory for lambda and outp */
         last_transition_matrix = allocateLambda(num_good_snps, 
                                         cliParams->num_haplotypes);
         last_emission_matrix = allocateOutp(num_good_snps, cliParams->num_haplotypes, 
                                     num_emission_types);
         readLamboutpFile(cliParams->lambda_outp_file, last_transition_matrix, 
                          last_emission_matrix, cliParams->num_haplotypes, 
                          num_emission_types, num_good_snps);
      }
      else
      {
         fprintf(stderr, "--sort_only or --no_hmm requires passing parameters with " 
                         "--lambda_outp_file\n");
         exit(EHMM_NOPARAM);
      }
      

      
      /* since the various path algorithms require us to pass a pointer to an 
         em_out_t structiure, we construct one here with necessary fields 
         initialized */
      emOutput = malloc(sizeof(em_out_t));
      emOutput->last_emission_matrix = last_emission_matrix;
      emOutput->last_transition_matrix = last_transition_matrix;
   
      emOutput->last_init_matrix = allocate2Dd(num_strains, 
                                              cliParams->num_haplotypes);   
      if (emOutput->last_init_matrix == NULL)
      { 
         fprintf(stderr, "Unable to allocate Init\n");
         exit(EHMM_NOMEM);
      } 
   
      for (int i = 0; i < num_strains; i++)
         for (int j = 0; j < cliParams->num_haplotypes; j++)
         {
            emOutput->last_init_matrix[i][j] = 1.0 / cliParams->num_haplotypes;
         }
         
         
      emOutput->last_log_liklihood = 0.0;
      dPredMatrix = allocatePred(num_good_snps + 2, cliParams->num_haplotypes);
      dFiltMatrix = allocatePred(num_good_snps + 2, cliParams->num_haplotypes);
      dSmthMatrix = allocatePred(num_good_snps + 2, cliParams->num_haplotypes);
      
      for (int i = 0; i < num_strains; i++)
      {
         filt(dPredMatrix, dFiltMatrix, SNPs, i, last_emission_matrix, 
              last_transition_matrix, emOutput->last_init_matrix, num_good_snps, 
              cliParams->num_haplotypes, cliParams->miss_option);
         emOutput->last_log_liklihood += loglike(SNPs, i, num_good_snps, 
                          last_emission_matrix, cliParams->num_haplotypes, 
                          dPredMatrix, num_emission_types);
      }   
       
      freePred(dPredMatrix, num_good_snps + 2);
      freePred(dFiltMatrix, num_good_snps + 2);
      freePred(dSmthMatrix, num_good_snps + 2);
      
      
      /* if user passed --no_hmm AND --graphviz, then we need to do a little 
      extra work. We need to compute the marginal state probability which is 
      used by the writeGraph function*/
      if (cliParams->graphviz)
      {
         marg_prob_matrix = allocate2Dd(num_good_snps, cliParams->num_haplotypes);

         if (marg_prob_matrix == NULL)
         {
            fprintf(stderr, 
                 "    Unable to allocate memory for Marginal Probability Matrix\n");
            exit(EHMM_NOMEM);
         }  
   
          computeMargProb(marg_prob_matrix, last_transition_matrix, 
                          cliParams->num_haplotypes, num_good_snps);
      
          printf("Writing graphviz file...");
        
          snprintf(path, MAXPATHLEN-1, "%sgraph.dot", file_prefix);
          writeGraph(last_transition_matrix, last_emission_matrix, marg_prob_matrix, num_good_snps,
                     cliParams->num_haplotypes, SNPs, 
                     num_emission_types, path, 100, 5);
          printf("done\n");
       
          free2Dd(marg_prob_matrix, num_good_snps);
      }
      
   }
     

   /* now we use the results to calculate a path and fill in the data */



   fill_confidences = allocate2Dd(num_good_snps, num_strains);
   
                               

   /* compute the simple path if necessary */
   if (cliParams->path_option == PATH_MAXSMTH || 
       cliParams->path_option == PATH_BOTH)
   {
      fill_and_save(SNPs, &emParams, emOutput, PATH_MAXSMTH, 
                    cliParams->sort_option, cliParams->confidence_threshold,
                    num_emission_types, cliParams->chromosome, 
                    cliParams->smth_out, fill_confidences);
   }
   
   
   /* compute viterbi path if necessary (this is the most common path algorithm
      and is the default behavior */
   if (cliParams->path_option == PATH_VITERBI ||
       cliParams->path_option == PATH_BOTH)
   {
      fill_and_save(SNPs, &emParams, emOutput, PATH_VITERBI,
                    cliParams->sort_option, cliParams->confidence_threshold,
                    num_emission_types, cliParams->chromosome, 
                    cliParams->smth_out, fill_confidences);
   }

   /* State Sorting  */

   if (cliParams->sort_option == SORT_MAX_TRACE)
   {

      sort_states(SNPs, num_good_snps, &emParams, emOutput,
                  last_transition_matrix, last_emission_matrix, 
                  cliParams->smth_out);  
      
   }
   
   temp_file_cleanup();

   free2Dd(fill_confidences, num_good_snps);
   freeLambda(last_transition_matrix, num_good_snps, cliParams->num_haplotypes);
   freeOutp(last_emission_matrix, num_good_snps, cliParams->num_haplotypes);
   freeEM_Out(emOutput, num_strains);
   
   /* if the user passed a prefix then we don't want to free file_prefix because
      it just points to an argv string */
   if (cliParams->user_prefix_flag == FALSE) {
      free(file_prefix);
   }
   
   fclose(fsSummaryFile);
  
   return 0;

}
