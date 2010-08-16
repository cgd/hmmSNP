/*
   File Name: dataFilter.c
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
      Algorithm designed by Jin Szatkeiwicz under the supervision of 
      Gary Churchill; and a matlab version of the algorithm is available from 
      Jin.
   Date developed: July-October 2006
   Purpose: This file contains the definition of the filterData function, 
     which is used by the program to remove "bad" SNPs
   Overview: filterData is called by the main program. This function removes 
     any SNP with a missingness that exceedes the maxmiss parameter, or that 
     fails the checks performed by the checkSNP function. checkSNP is defined 
     locally to this file.
   Usage: filterData function is described bellow
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

#include "constants.h"
#include "dataIO.h"





/*******************************************************************************
 * filterData  creates list of "good SNPs" so the nth snp_t struct in the list
 *   is the nth good SNP not the nth overall SNP
 *
 * INPUTS:  snp_t **filtered_snps: array of pointers to snp_t structs
 *          snp_t *snps: original unfiltered data (array of snp_t structs)
 *          int num_snps: number of SNPs in original data set
 * RETURNS: void
 *          returns pointers to good SNPs via  filtered_snps
 *
 * ASSUMES: **filtered_snps has been allocated
 *
 * EFFECTS: alters contents of filtered_snps
 *
 * ERROR CONDITIONS: none
 *
 * COMMENTS: filtered_snps elements are pointers into rows into snps, 
 *    so we don't actually copy data.  if you iterate through filtered_snps
 *    you will skip over rows in snps that have missing > dMaxMiss or 
 *    other problems
 *
 ******************************************************************************/
void filterData(snp_t **filtered_snps, snp_t *snps, int num_snps)
{

   int good_counter = 0;
   
   for (int i = 0; i < num_snps; i++)
   {
      /* if there is a problem with the SNP, other than missingness, 
         flag it by setting the missingness to 100% */
         
      if (snps[i].iStatus == SNP_GOOD)
      {
         filtered_snps[good_counter++] = &snps[i];
      }
      else 
      {
         /* set haplotypes to zero to flag the unfiltered data so we know 
            this SNP has been removed*/
         snps[i].iNumHaplotypes = 0;
      }
   }

}

/*******************************************************************************
 * double getMissingness(snp_t *psnp) - calculates missing % of SNP
 *
 * INPUTS: pointer to a snp_t struct
 *
 * RETURNS: double precisions missingness
 *
 * ASSUMES: psnp points to initialized snp_t struct
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
double getMissingness(snp_t *psnp)
{
   int num_strains;
   int num_missing = 0;
   
   num_strains = getNumStrains();
   
   for (int i = 0; i < num_strains; i++)
   {
      if (psnp->iEncodedSNP[i] == 2)
      {
         num_missing++;
      }
   }
   
   return (double)num_missing / (double)num_strains;
}
