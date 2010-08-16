/*
   File Name: graph.c
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
      Algorithm designed by Jin Szatkeiwicz under the supervision of 
      Gary Churchill; and a matlab version of the algorithm is available from 
      Jin.
   Date developed: July-October 2006
   Purpose: contains functions related to creating graphviz output
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

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>


#include "constants.h"
#include "graph.h"
#include "dataIO.h"

int *nodeColors(double**, int, int);



/*******************************************************************************
 * writeGraph - writes out graphviz file
 *
 * INPUTS:  double ***transition_matrix - transition matrix used to draw edges
 *          double ***emission_matrix - emission matrix, used for node labels
 *          int num_snps - number of SNPs
 *          int num_haplotypes - number of different haplotypes
 *          snp_t *snps - array of snp_t structs. used for Allele 
 *            for each SNP
 *          int miss - how many possible output options (2 or 3) depending on 
 *            missing as emission or random
 *          const char path[] - output file path
 *          int width - graphviz bounding box width (inches)
 *          int height - graphviz bounding box height (inches )
 *
 * RETURNS: void
 *
 * ASSUMES:
 * 
 * EFFECTS:
 *
 * ERROR CONDITIONS: will exit if unable to open file.
 *
 * COMMENTS:  num_haplotypes will have to be changed to an array if we allow 
 *            more than one haplotype per SNP
 *
 ******************************************************************************/
void writeGraph(double ***transition_matrix, double ***emission_matrix, 
                double **marg_prob_matrix, int num_snps,
                int num_haplotypes, snp_t *snps, int miss, 
                const char path[], int width, int height)
{
   FILE *graph_fs;
   int max;
   double probability;
   double line_weight;
   char label;
   int node_id;
   int target_node;
   int *color_array;
   
   if (HMM_VERBOSITY > 0) printf("attempting to open %s for writing\n", path);
   
   graph_fs = fopen(path, "w");
   
   if (graph_fs == NULL)
   {
      perror("Unable to open graphviz output file");
      exit(errno);
   }
   
   color_array = nodeColors(marg_prob_matrix, num_haplotypes, num_snps);

   
   /* print graphviz header */
   fprintf(graph_fs, "digraph G {\n");
   fprintf(graph_fs, "center = 1;\n");
   fprintf(graph_fs, "size=\"%d,%d\";\n", width, height);
   fprintf(graph_fs, "rankdir=LR;\n");
   fprintf(graph_fs, "node [shape=circle];\n");
   fprintf(graph_fs, "edge [arrowhead=none];\n");
   
   
   /* print nodes, start with the "begin state" */
   fprintf(graph_fs, "1 [ label=\"b\" style=\"filled,setlinewidth(3.0)\" fillcolor=grey90 ];\n");
   node_id = 2;
   for (int i = 0; i < num_snps; i++)
   {  
      /* find most likely output for each hap, and use for node label/color */
      for (int j = 0; j < num_haplotypes; j++)
      {
         max = 0;
         for (int k = 1; k < miss; k++)
         {
            if (emission_matrix[i + 1][j][k] > emission_matrix[i + 1][j][max])
            {
               max = k;
            }
         }
 
         probability = emission_matrix[i + 1][j][max];
         
         if (max < 2)
         {
            label = snps[i].cAlleleArray[max];
         }
         else
         {
            label = 'N';
         }
          
         if (marg_prob_matrix[i][j] != 0.0)
         {
            fprintf(graph_fs,
                    "%d [ label=\"%c\" style=\"filled,setlinewidth(3.0)\" fillcolor=grey%d ];\n",
                    node_id, label, color_array[node_id-2]);
         }
         else
         {
         /*
            fprintf(graph_fs,
                    "%d [ style=filled color=grey%d ];\n",
                    node_id, color_array[node_id-2]);
                    */
         }
         ++node_id;
         
         
      }
   }
   fprintf(graph_fs,
           "%d [ label=\"e\" style=\"filled,setlinewidth(3.0)\" fillcolor=grey90 ];\n", node_id);

   

   /* print edges */

   node_id = 2;
   target_node = 2;
   
   /* from b */
   for (int i = 0; i < num_haplotypes; i++)
   {
      fprintf(graph_fs, "1 -> %d [style=\"setlinewidth(%2.1f)\"];\n", target_node++, 
              transition_matrix[0][0][i+1]*5);
   }

   
   /* from other nodes */
   for (int i = 1; i < num_snps; i++)
   {
      for (int j = 0; j < num_haplotypes; j++)
      {
         for (int k = 0; k < num_haplotypes; k++)
         {
            line_weight = transition_matrix[i][j][k]*5;
            if (line_weight >= EDGE_THRESHOLD)
            {
               fprintf(graph_fs, 
                       "%d -> %d [style=\"setlinewidth(%2.1f)\"];\n", 
                       node_id, target_node, line_weight);
            }
            else if (transition_matrix[i][j][k] != 0.0)
            {
               fprintf(graph_fs, "%d -> %d [style=invis];\n", node_id, 
                       target_node);            
            }
            ++target_node;
         }
         target_node -= num_haplotypes;
         ++node_id;
      }
      target_node += num_haplotypes;
   }
   
   /* to e */
   
   for (int i = 0; i < num_haplotypes; i++)
   {
      if (transition_matrix[num_snps][i][0] > 0.0)
      {
         fprintf(graph_fs, "%d -> %d [style=\"setlinewidth(%2.1f)\"];\n", 
                 node_id, target_node, transition_matrix[num_snps][i][0]*5);  
      }
      ++node_id;
   }
   
   fprintf(graph_fs, "}\n");
   
   fclose (graph_fs);
   
}

/*******************************************************************************
 * nodeColors - compute and return node colors
 *     
 * INPUTS: double **marg_prob_matrix - marginal state probabilities
 *         int num_haplotypes - number of different haplotypes
 *         int num_snps - number of SNPs
 *
 * RETURNS: pointer to integer array of color intensities on a scale of 0 to 100
 *
 * ASSUMES:
 *
 * EFFECTS: allocates memory for node color intensities
 *
 * ERROR CONDITIONS: could cause segfault if for some reason malloc failed 
 *   (should trap and exit cleanly)
 *
 * COMMENTS:
 *
 ******************************************************************************/
int *nodeColors(double **marg_prob_matrix, int num_haplotypes, int num_snps)
{
   int node_id;
   int *color_array;

   color_array = (int*)malloc((num_haplotypes * num_snps) * sizeof(int)); 

   
   node_id = 0;
   for (int t = 0; t < num_snps; t++)
   {
      for (int i = 0; i < num_haplotypes; i++)
      {
         color_array[node_id++] = (int)((1 - marg_prob_matrix[t][i]) * 100.0 + 0.5);

      }
   }
   return color_array;
}
