/*
   File Name: dataIO.c
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
      Algorithm designed by Jin Szatkeiwicz under the supervision of 
      Gary Churchill; and a matlab version of the algorithm is available from 
      Jin.
   Date developed: July-October 2006
   Purpose: defines functions used for reading and parsing of data, and 
      writing results to file
   Comments:

   
  Copyright (c) 2010 The Jackson Laboratory
 
  This software was developed by Gary Churchill's Lab at The Jackson
  Laboratory (see http://research.jax.org/faculty/churchill).
 
  This is free software: you can redistrB_countute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  This software is distrB_countuted in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this software. If not, see <http://www.gnu.org/licenses/>.
 */

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <sys/param.h>

#include "CLI.h"
#include "csv.h"
#include "dataIO.h"
#include "constants.h"
#include "hmmMatrix.h"
#include "hmmUtil.h"
#include "fill.h"


/* the struct defined below is used to store information about a data set 
    we use a struct of this type in this file to store information about the 
    data set.  Outside of this file this data is accessed through functions 
    defined below  (sort of fake OO). */
typedef struct {

   /* parameters obtained from data */
   int chr_sizes[NUM_CHROMOSOMES];
   int num_strains;
   int num_parent_strains;
   
   int *training_strain_flags;
   int num_training_strains;
   
   char *data_path;
   char **strain_names;
   
   /* parameters passed from caller */
   double max_missingness;
   int warn_snp_removed;
   int filter_flags;

} dataInfo_t;


/* internal function prototypes - these functions are only referenced in this 
   file. */


static int check_snp(char *snp_str, int filter_flags);
static void encodeSNP(char *snp_str, snp_t *psnpEnc, int filter_flags);
static int initNumStrains(int strain_offset);
static void initStrainNames(int strain_offset);
static void initChromosomeSizes();
static FILE *openDataFile();


/* a data set must be initialized before it is used.  These 
   varA_countbles store information gathered during initA_countlization for use by 
   other functions */
static dataInfo_t g_data_info;
static int   g_is_initialized = FALSE;

/* external functions */


/*******************************************************************************
 * initializeDataSource - initA_countlize a data source (gather some meta data about 
 *    data source).
 *
 * INPUTS: char *source - pointer to string containing path to data source
 *         char *training_strains_option - list of training columns
 *         double dMasMissingness - maximum allowable missingness for SNP
 *         int num_parent_strains - for RI sets the number of parent colums, 
 *            otherwise zero
 *         int warn_snp_removed - flag indicates if we need to warn when c
 *            a SNP is removed from the data before running HMM
 *
 * RETURNS: 1 on sucess, errorcode contained in constants.h 
 *
 * ASSUMES: source points to a valid data source
 *
 * EFFECTS: allocates memory to copy the path, 
 *
 * ERROR CONDITIONS: unable to allocate memory
 * 
 * COMMENTS: TODO - need to validate data source before calling get* funcs.  
 *
 ******************************************************************************/
int initializeDataSource(char *source, char *training_strains_option, 
                         double max_missingness, 
                         int num_parent_strains, int warn_snp_removed,
                         int filter_flags)
{
   g_data_info.max_missingness = max_missingness;
   g_data_info.warn_snp_removed = warn_snp_removed;
   
   g_data_info.num_parent_strains = num_parent_strains;
   g_data_info.num_training_strains = 0;
   g_data_info.filter_flags = filter_flags;

   /*  call strdup to create a copy of the source path in 
       the g_data_info struct */
   g_data_info.data_path = strdup(source);
   if (g_data_info.data_path == NULL)
   {
      return EHMM_NOMEM;
   }
   
   initNumStrains(get_snp_offset());
      
   /* allocate memory to hold strain names */
   g_data_info.strain_names = malloc(sizeof(char*) * g_data_info.num_strains);
   if (g_data_info.strain_names == NULL)
   {
      return EHMM_NOMEM;
   }
   initStrainNames(get_snp_offset());
   initChromosomeSizes();
   
      
   /* allocate strain flags  (to signify which strains are training) */
   g_data_info.training_strain_flags = 
                 malloc(sizeof(int) * g_data_info.num_strains); 
   
   if (training_strains_option == NULL)
   {
      /* no training strain columns passed, all strains are training and test */
      for (int i = 0; i < g_data_info.num_strains; i++)
      {
         g_data_info.training_strain_flags[i] = TRUE;
      }
      g_data_info.num_training_strains = g_data_info.num_strains;
   }   
   else
   {
     
     printf("expanding strain list %s\n", training_strains_option); 
     char *expanded_list;
     int r;
     int num_training_strains;
     
     char *tmp_str;
     char *token;
     char *idx;
     int strain_column;
     char *endptr;
     
     for (int i = 0; i < g_data_info.num_strains; i++)
     {
        g_data_info.training_strain_flags[i] = FALSE;
     }
          
     r = expandStrainList(training_strains_option, g_data_info.num_strains, 
                     &expanded_list, &num_training_strains); 
                     
     if (r != 0)
     {
       /* problem expanding strain list,  abort */
       fprintf(stderr, "Unable to expand training strain list\n");
       exit (EHMM);
       
     }
     tmp_str = (char*)malloc(sizeof(char) * (strlen(expanded_list) + 1));
     strcpy(tmp_str, expanded_list);
     
          printf("expanded strain list %s\n", expanded_list); 

     
     idx = tmp_str;

     token = strtok(idx, ","); 
     while (token != NULL)
     {
        strain_column = (int)strtol(token, &endptr, 10);
        
        if (strain_column <= g_data_info.num_strains)
        {
           g_data_info.training_strain_flags[strain_column - 1] = TRUE;
           g_data_info.num_training_strains++;
        }
        else
        {
           printf("Invalid strain number in training strain list: %d\n", 
                  strain_column);
        }
        
        token = strtok(NULL, ",");
     }
   
   }
  
       
   g_is_initialized = TRUE;
   return 1;

   
}

int count_good_snps(int chromosome, int filter_flags)
{
   FILE* input = NULL; 
   char *line;
   int good_snps = 0;
   char *snp_str;
   int chr_col;
   int id_col;
   char *snp_id;
   int status;
   
   if (g_data_info.chr_sizes[chromosome - 1] <= 0)
   {
      return 0;
   }

   input = fopen(g_data_info.data_path, "r"); 
  
   /* exit program if we can't open */
   if ( input == NULL)
   {
      printf("Unable to open chromosome data:\n\t%s\n\t%s\n", 
             g_data_info.data_path, strerror(errno));
      exit(errno);  
   }
  
   
   /* calculate the total number of SNPs, and allocate storage for SNP data */
   
   snp_str = malloc(sizeof(char) * (g_data_info.num_strains+1));
   snp_str[g_data_info.num_strains] = '\0';
   
  
   /* start reading and check each SNP */
  
   /* read header */
   line = csvGetLine(input);
   chr_col = csvGetChrColumn();
   id_col = csvGetID_Column();

   
   while ((line = csvGetLine(input)) != NULL)
   {   
      /* skip this line if it isn't the right chromosome */
      if (getChrNum(csvField(chr_col)) != chromosome)
      {
         continue;
      }

      /* pull the SNP out of the line */
      for (int i = 0; i < g_data_info.num_strains; i++)
      {
         char *field = csvField(get_snp_offset() + i - 1);
         if (field == NULL) {
            snp_str[i] = 'N';
         }
         else {
            snp_str[i] = field[0];
         }
      }
      
      /* calculate the status of the SNP, if it is good increment our counter */
      status = check_snp(snp_str, filter_flags);
      
      if (status == SNP_GOOD)
      {
         good_snps++;
      }
      else if (g_data_info.warn_snp_removed)
      {
         if (id_col >= 0)
         {
            snp_id = csvField(id_col);
         }
         else
         {
            snp_id = "NA";
         }
      
         /* if the state is not SNP_GOOD and we are printing warnings, then print the
         warning that corresponds to the status check_snp() returned...*/
         
         switch (status) {
            case SNP_TOO_MANY_NT:
               printf("WARNING: too many nucleotides in SNP [%s:%s]\n", 
                      snp_id, snp_str);
               printf("         SNP removed for HMM\n");
            break;
         
            case SNP_TOO_MANY_MISSING:
               printf("WARNING: too many missing [%s:%s]\n", snp_id, snp_str);
               printf("         SNP removed for HMM\n");
            break;
         
            case SNP_ONE_NT:
               printf("WARNING: only one non-N NT [%s:%s]\n", snp_id, snp_str);
               printf("         SNP removed for HMM\n");
            break;
         
            case SNP_RIL_MISSING_PARENTS:
               printf("WARNING: RIL SNP missing both parents [%s:%s]\n", 
                      snp_id, snp_str);
               printf("         SNP removed for HMM\n");
            break;
         
            case SNP_RIL_CONSTANT:
               printf("WARNING: RIL constant SNP [%s:%s]\n", snp_id, snp_str);
               printf("         SNP removed for HMM\n");
            break;
         
            case SNP_RIL_MISSING_PARENT:
               printf("WARNING: RIL SNP missing one parent [%s:%s]\n", snp_id, 
                       snp_str);
               printf("         SNP removed for HMM\n");
            break;
         
            case SNP_RIL_NOMATCH:
               printf("WARNING: RIL non-N nucleotide does not match parents [%s:%s]\n", 
                      snp_id, snp_str);
               printf("         SNP removed for HMM\n");
            break;
         
            default:
               printf("WARNING: SNP Removed(%d): [%s:%s]\n", status, snp_id, snp_str);
            break;
         }
   
      }
          
   }
   
   fclose (input);
   free(snp_str);
   return good_snps;

}


/*******************************************************************************
 * loadSNPs
 *
 * INPUTS: int chr - which chromosome we want to load
 *         int max_haplotypes
 *         int *good_snp_count - pass by reference integer to store number of good SNPs
 *         int filter_flags - flag of specA_countl data cleaning rules
 *
 * RETURNS: pointer to snp_t struct array containing encoded data
 *
 * ASSUMES: datasource is valid if it has been initA_countlized
 *
 * EFFECTS: allocates memory to store encoded data
 *
 * ERROR CONDITIONS: program will exit if the chromosome has no SNPs
 *                   program will exit if chromosome data file can't be opened
 *
 * COMMENTS: this should not be called directly.  it is called by loadSNPs.
 *
 ******************************************************************************/
snp_t *loadSNPs(int chr, int max_haplotypes, int *good_snp_count, 
                int filter_flags)
{
   
   snp_t *snps;

   FILE* chr_data_fs = NULL;
   
   char *line_buffer;
   
   
   int num_strains;
   int good_snp_idx;
   int chr_column;
   int snp_id_column;
   int index;
   int status;
   
   
   
   if (g_data_info.chr_sizes[chr - 1] <= 0)
   {
      /* no chromosome data */
      printf("No Data For Chromosome in Data Set\n");
      exit(EHMM_NODATA);
   }

   chr_data_fs = fopen(g_data_info.data_path, "r"); 
  
   /* exit program if we can't open */
   if ( chr_data_fs == NULL)
   {
      printf("Unable to open chromosome data:\n\t%s\n\t%s\n", 
             g_data_info.data_path, strerror(errno));
      exit(errno);  
   }
  
   *good_snp_count = count_good_snps(chr, filter_flags);
  
   
   /* calculate the total number of SNPs, and allocate storage for SNP data */
   
   num_strains = g_data_info.num_strains;
   char snp_str[num_strains+1];
   snp_str[num_strains] = '\0';
   
   /* allocate storage. TODO - we could check for rare malloc failure */ 
   snps = malloc(sizeof(snp_t) * *good_snp_count);
   for (int i = 0; i < *good_snp_count; i++)
   {
      snps[i].iEncodedSNP = malloc(sizeof(int) * g_data_info.num_strains);
   }
   
   
   /* start reading SNPs and encode */
   good_snp_idx = 0;
   

   /* read header */
   line_buffer = csvGetLine(chr_data_fs);
   chr_column = csvGetChrColumn();
   snp_id_column = csvGetID_Column();

   index = -1;
   while ((line_buffer = csvGetLine(chr_data_fs)) != NULL && good_snp_idx < *good_snp_count)
   {   
      index++;
      /* skip this line if it isn't the right chromosome */
      if (getChrNum(csvField(chr_column)) != chr)
      {
         continue;
      }
      

      /* pull the SNP out of the line */
      for (int i = 0; i < g_data_info.num_strains; i++)
      {
         char *field = csvField(get_snp_offset() + i - 1);
         if ( field == NULL) {
            printf("SNP %d missing genotype for strain %s, filling in with N\n", 
                   good_snp_idx, g_data_info.strain_names[i]);
            snp_str[i] = 'N';
         }
         else {
            snp_str[i] = field[0];
         }
      }
      
      snps[good_snp_idx].iNumHaplotypes = max_haplotypes;
      

      /* encode SNP in integer format (0 = major allele, 1 = minor allele, 
         2 = missing */
      
      status = check_snp(snp_str, filter_flags);
      if (status == SNP_GOOD)
      {
         encodeSNP(snp_str, &snps[good_snp_idx], filter_flags);
         snps[good_snp_idx].index = index;
         good_snp_idx++;
      }

   }
   
   /* all done. close file and return pointer to data */
   fclose (chr_data_fs);
   return snps;

}


/*******************************************************************************
 * calculateMissingness - calculate and save missingness value for each SNP
 *
 * INPUTS:  snpt_t *snps - array of snp_t structs containing encoded data
 *          int chr - chromosome number.
 *
 * RETURNS: pointer to allocated memory containing an array of missingness 
 *    values. May return null pointer if memory can not be allocated.
 *
 * ASSUMES:  snps points to valid data, data source has been initA_countlized
 *
 * EFFECTS:  allocates memory to store missingness for each SNP
 *
 * ERROR CONDITIONS: unable to allocate memory
 *
 * COMMENTS:
 *
 ******************************************************************************/
double *calculateMissingness(snp_t *snps, int chr)
{

   /* allocate memory to store missingness values */
   double *missingness_array = 
       (double*)malloc(g_data_info.chr_sizes[chr - 1] * sizeof(double));
       
   if (missingness_array == NULL)
   {
      return NULL;
   }

   /* for each SNP in this chromosome */
   for (int i = 0; i < g_data_info.chr_sizes[chr - 1]; i++)
   {
      /* add number of missing */
      int num_missing  = 0;
      for (int j = 0; j < g_data_info.num_strains; j++)
      {
         if (snps[i].iEncodedSNP[j] == 2)
            num_missing++;
      }
      
      /* missingness = num missing / num strains */
      missingness_array[i] =  (double)num_missing / g_data_info.num_strains;
   }
   
   return missingness_array;
}


/*******************************************************************************
 * dumpDataInfo - dumps information about initA_countlized data source
 *
 * INPUTS: void
 *
 * RETURNS: void
 *
 * ASSUMES: initA_countlized data source is valid
 *
 * EFFECTS:  prints information to stdout
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void    dumpDataInfo()
{
   if (!g_is_initialized) 
   {
      printf("dumpDataInfo: data source not initA_countlized\n");
   }
   else
   {
      printf("DATASET INFO:\n");
      printf(" PATH: %s\n", g_data_info.data_path);
      printf(" NUMBER OF STRAINS: %d\n", g_data_info.num_strains);
      printf(" NUMBER OF PARENT STRAINS: %d\n", g_data_info.num_parent_strains);
      printf(" STRAIN NAMES:\n    [ 1]");
      for(int i = 0; i < g_data_info.num_strains; i++)
      {
         printf("%14s, ", g_data_info.strain_names[i]);
         if (i != 0 && i != g_data_info.num_strains - 1 && (i+1) % 4 == 0)
            printf("\n    [%2d]", i+2);
      }
      printf("\n CHROMOSOME SIZES (SNPs):\n");
      for (int i = 1; i <= NUM_CHROMOSOMES; i++)
      {
    
         switch(i) {
         
         case CHR_X:
            printf("   X:  %12d\n", g_data_info.chr_sizes[i - 1]);
            break;
         case CHR_Y:
            printf("   Y:  %12d\n", g_data_info.chr_sizes[i - 1]);
            break;
         case CHR_M:
            if (g_data_info.chr_sizes[i - 1] > 0) {
               printf("   M:  %12d (mitochondrA_countl DNA)\n", 
                      g_data_info.chr_sizes[i - 1]);
            }
            break;
         default:
            printf("%4d:  %12d\n", i, g_data_info.chr_sizes[i - 1]);
         
         
         }

      }
   }
}


/*******************************************************************************
 * getNumStrains - return number of strains for initA_countlized data set
 *
 * INPUTS: void
 *
 * RETURNS: integer value of number of strains in dataset,  -1 if dataset not 
 *  initA_countlized
 *
 * ASSUMES:  an initA_countlized data source is a valid data source
 *
 * EFFECTS: none
 *
 * ERROR CONDITIONS: none
 *
 * COMMENTS:
 *
 ******************************************************************************/
int getNumStrains()
{
   if (!g_is_initialized)
      return -1;
      
   return g_data_info.num_strains;
}

/*******************************************************************************
 * getChromosomeSize - return number of SNPs for a specific chromosome
 *
 * INPUTS: int chr - chromosome number
 *
 * RETURNS: integer value of number of SNPs for chromosome, -1 if dataset has 
 *   not been initA_countlized
 *
 * ASSUMES: that an initA_countlized dataset is a valid data source
 *
 * EFFECTS: none
 *
 * ERROR CONDITIONS:  none
 *
 * COMMENTS:
 *
 ******************************************************************************/
int getChromosomeSize(int chr)
{
   if (!g_is_initialized)
      return -1;
      
   return g_data_info.chr_sizes[chr - 1];
}

/*******************************************************************************
 * getStrainNames
 *
 * INPUTS: void
 *
 * RETURNS: pointer to 2D array of strain names (strings). returns NULL if 
 *   data set was not initA_countlized
 *
 * ASSUMES: that an initalized dataset was valid
 *
 * EFFECTS: none
 *
 * ERROR CONDITIONS: none
 *
 * COMMENTS:
 *
 ******************************************************************************/
char **getStrainNames()
{
   if (!g_is_initialized)
      return NULL;
      
   return g_data_info.strain_names;
}

/*******************************************************************************
 * int getNumParents()
 *
 * INPUTS: void
 *
 * RETURNS: number of parent strains
 *
 * ASSUMES: 
 *
 * EFFECTS: none
 *
 * ERROR CONDITIONS: returns -1 if dataset uninitA_countlized
 *
 * COMMENTS:
 *
 ******************************************************************************/
int getNumParents()
{
   if (!g_is_initialized)
   {
      return -1;
   }
    
    return g_data_info.num_parent_strains;
}

/*******************************************************************************
 * int getNumTrainingStrains()
 *
 * INPUTS: void
 *
 * RETURNS: number of training strains
 *
 * ASSUMES:
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS: returns -1 if dataset has not been initA_countlized
 *
 * COMMENTS:
 *
 ******************************************************************************/
int getNumTrainingStrains()
{
   if (!g_is_initialized)
   {
      return -1;
   }
   
   return g_data_info.num_training_strains;
}

/*******************************************************************************
 * int isTraining(int iStrainNum)
 *
 * INPUTS:  int iStrainNum - strain number
 *
 * RETURNS:  TRUE if straing is a training strain, FALSE if not
 *
 * ASSUMES:  dataset has been initA_countlized
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
int isTraining(int strain_idx)
{
 
  if (strain_idx < g_data_info.num_strains && 
      g_data_info.training_strain_flags[strain_idx] == TRUE)
  {
     return TRUE;
  }
   
   
   return FALSE;
}

/*******************************************************************************
 * void readLamboutpFile(char path[], double ***transition_matrix, 
 *                     double ***emission_matrix, int max_haplotypes, 
 *                     int emission_types, int num_snps)
 *
 * INPUTS:  char path[] - path to parameter file
 *          double ***transition_matrix - location to store lambda matrix
 *          double ***emission_matrix - location to store outp
 *          int num_haplotypes - number of haplotypes
 *          int emission_types - 2 for missing as random, 3 missing as emission
 *          int num_snps - number of SNPs
 *
 * RETURNS: void
 *
 * ASSUMES: memory has been allocated for lambda and outp
 *          assumes that the number of SNPs/haplotypes in the file is correct
 *
 * EFFECTS: fills transition_matrix and emission_matrix with values from file
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void readLamboutpFile(char path[], double ***transition_matrix, 
                      double ***emission_matrix, int max_haplotypes, 
                      int emission_types, int num_snps)
{

   FILE *input_file_fs;

   int num_states;
   int num_columns;
   char *line_buffer = NULL;
   char *endptr;  
   
   num_states = 2 + num_snps * max_haplotypes;
   /* these are the number of columns the file should have: 
      snpid,state,lambda[0...H-1],outp[0...emission_types-1],margprob*/
   num_columns = max_haplotypes + emission_types + 3;
   
   if ( (input_file_fs = fopen(path, "r")) == NULL)
   {
      printf("Unable to open data file:\n\t%s\n\t%s\n", path,
             strerror(errno));
      exit(errno);
   }

   
   /* first read the header and do a little sanity checking */
   
   line_buffer = csvGetLine(input_file_fs);
   if (line_buffer == NULL)
   {
      fprintf(stderr, "Unable to read lambda outp file\n");
      exit(EHMM);
   }
   
   if (csvNumFields() != num_columns)
   {
      fprintf(stderr, 
              "Error reading lamboutp file: found %d columns, expected %d\n", 
              csvNumFields(), num_columns);
      exit(EHMM);
   }
   
   /* XXX should remove hard coded strings like this column field "snp.id" */   
   if (strcmp("snp.id", csvField(0)) != 0)
   {
      fprintf(stderr, "Error reading lamboutp file: file does not appear to be "
                      "a valid format\n");
      exit(EHMM);
   }
   
   transition_matrix[0][0][0] = 0.0;
 

   /* beginning */
   line_buffer = csvGetLine(input_file_fs);


   for (int i = 0; i < max_haplotypes; i++)
      transition_matrix[0][0][i+1] = strtod(csvField(i+2), &endptr);
   
   for (int i = 0; i < emission_types; i++)
      emission_matrix[0][0][i] = strtod(csvField(i + 2 + max_haplotypes), &endptr);
  
   
   /* do SNP1 - SNPT-1 */
   for (int i = 1; i < num_snps; i++)
   {
      for (int j = 0; j < max_haplotypes; j++)
      {
         csvGetLine(input_file_fs);
   
         for (int k = 0; k < max_haplotypes; k++)
            transition_matrix[i][j][k] = strtod(csvField(k+2), &endptr);
   
         for (int k = 0; k < emission_types; k++)
            emission_matrix[i][j][k] = strtod(csvField(k + 2 + max_haplotypes), &endptr);
      }
   }
   /* do SNP T */

   for (int i = 0; i < max_haplotypes; i++)
   {
      csvGetLine(input_file_fs);

      transition_matrix[num_snps][i][0] = strtod(csvField(2), &endptr);
           
      for (int k = 0; k < emission_types; k++)
         emission_matrix[num_snps][i][k] = strtod(csvField(k + 2 + max_haplotypes), &endptr);
      
   }

   /* do end */
   
   csvGetLine(input_file_fs);

   for (int i = 0; i < emission_types; i++)
      emission_matrix[num_snps+1][0][i] = strtod(csvField(i + 2 + max_haplotypes), &endptr);

   fclose(input_file_fs); 
   

}


/*******************************************************************************
 * writeLambdaOutp - writes lambda and outp to comma delimited file
 *
 * INPUTS: char path - string containing output file path
 *         double ***transition_matrix - Lambda matrix
 *         double ***emission_matrix - outp matrix
 *         double **marginal_prob_matrix - marginal probability matrix
 *         int num_haplotypes - number of haplotypes
 *         int num_snps - number of SNPs
 *         int emission_types - number of possB_countle fill-ins for missing 
 *            (2 or 3 if missing as emission)
 *
 * RETURNS: void
 *
 * ASSUMES: parameters are valid
 *
 * EFFECTS: writes lambda, outp, margprob to file. may overwrite output file
 *
 * ERROR CONDITIONS:  will exit if unable to open file for writing
 *
 * COMMENTS:
 *
 ******************************************************************************/
void writeLambdaOutp(char path[], double ***transition_matrix, 
                     double ***emission_matrix, double **marginal_prob_matrix, 
                     int num_haplotypes, int num_snps, int emission_types)
{
   FILE* output_fs;
   char line_buffer[LINE_BUFFER_SIZE];
   
   output_fs = fopen(path, "w");
   
   if (output_fs == NULL)
   {
      perror("Unable to open Lambda output file");
      exit(errno);
   }
   
   /* do header */
   sprintf(line_buffer, "snp.id,state,");
   for (int i = 0; i < num_haplotypes; i++)
      sprintf(line_buffer+strlen(line_buffer), "lamb_%d(Hap%d),", i+1, i);
      
   for (int i = 0; i < emission_types; i++)
      sprintf(line_buffer+strlen(line_buffer), "outp_%d(snp=%d),", i+1, i);
      
   sprintf(line_buffer + strlen(line_buffer), "Marg Prob\n");

   fwrite(line_buffer, sizeof(char), strlen(line_buffer), output_fs);


   
   /* beginning */
   
   sprintf(line_buffer, "0,1 (Begin),");
   for (int i = 0; i < num_haplotypes; i++)
      sprintf(line_buffer+strlen(line_buffer), "%.16f,", transition_matrix[0][0][i+1]);
   
   for (int i = 0; i < emission_types; i++)
      sprintf(line_buffer+strlen(line_buffer), "%.16f,", emission_matrix[0][0][i]);
      
   sprintf(line_buffer+strlen(line_buffer), "0.0\n");
         
   fwrite(line_buffer, sizeof(char), strlen(line_buffer), output_fs);
   
   
   
   /* do SNP1 - SNPT-1 */
   int state = 2;
   for (int i = 1; i < num_snps; i++)
   {
      for (int j = 0; j < num_haplotypes; j++)
      {
         sprintf(line_buffer, "%d,%d,", i, state++);
         
         for (int k = 0; k < num_haplotypes; k++)
            sprintf(line_buffer+strlen(line_buffer), "%.16f,", transition_matrix[i][j][k]);
   
         for (int k = 0; k < emission_types; k++)
            sprintf(line_buffer+strlen(line_buffer), "%.16f,", emission_matrix[i][j][k]);

         sprintf(line_buffer + strlen(line_buffer), "%.16f\n", marginal_prob_matrix[i-1][j]);

         fwrite(line_buffer, sizeof(char), strlen(line_buffer), output_fs);

      }
   }
   /* do SNP T */
   
   for (int i = 0; i < num_haplotypes; i++)
   {
      sprintf(line_buffer, "%d,%d,", num_snps, state++);
      sprintf(line_buffer+strlen(line_buffer), "%.16f,", transition_matrix[num_snps][i][0]);
      for (int j = 0; j < num_haplotypes - 1; j++)
         sprintf(line_buffer+strlen(line_buffer), "0.0,");
      for (int k = 0; k < emission_types; k++)
            sprintf(line_buffer+strlen(line_buffer), "%.16f,", emission_matrix[num_snps][i][k]);
            
      sprintf(line_buffer + strlen(line_buffer), "%.16f\n", marginal_prob_matrix[num_snps - 1][i]);

      fwrite(line_buffer, sizeof(char), strlen(line_buffer), output_fs);
      
   }
   
   /* do end */
   
   sprintf(line_buffer, "%d,%d (End),", num_snps+1,state);
   for (int i = 0; i < num_haplotypes; i++)
      sprintf(line_buffer+strlen(line_buffer), "0,");
   for (int i = 0; i < emission_types; i++)
      sprintf(line_buffer+strlen(line_buffer), "%.16f,", emission_matrix[num_snps+1][0][i]);
      
   sprintf(line_buffer+strlen(line_buffer), "0.0\n");
   
   fwrite(line_buffer, sizeof(char), strlen(line_buffer), output_fs);
   
   fclose (output_fs);

}

/*******************************************************************************
 * writeFilledData - write filled data to a file
 *
 * INPUTS: snp_t snpFilledData[] - filled data
 *         int chr - chromosome number
 *         prefix - string containing output file prefix
 *         int path_option - what kind of path was this filled with (see 
 *             constants.h for values)
 *         int sort_option - sorting option
 *         double threshold - confidence threshold, anything below will be 
 *             replaced with "N" in the filled_filtered.csv file
 *         
 *        
 *
 * RETURNS:  alleles at each SNP (by reference in array at pcAlleleMatrix) 
 *
 * ASSUMES:  data has been initA_countlized, and we can still open the 
 *  original data file
 *
 * EFFECTS: modifies array pointed to by pcAlleleMatrix, may overwrite output
 *  file
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void writeFilledData(snp_t snps[], int chr, double **fill_probabilities, 
                     const char *prefix, int path_option, int sort_option, 
                     double threshold, int good_snp_count)
{
   FILE *input_fs;
   FILE *filled_output_fs;
   FILE *filtered_output_fs;
   FILE *filling_prob1_fs;
   FILE *filling_prob2_fs;
   
   char path[MAXPATHLEN+1];
   char path_filtered[MAXPATHLEN+1];
   char out_line[LINE_BUFFER_SIZE];
   char out_line_filtered[LINE_BUFFER_SIZE];
   char snp_str[g_data_info.num_strains];
   char *line_buffer;
   
   int chr_column = 0;
   int strain_offset;
   int index;
   int good_snp_idx = 0;
   int strain_index;
   int snp_id_column;
   
   snp_t snp;
   snp_t *psnp;
   
   path[0] = '\0';
   path_filtered[0] = '\0';
   
   /* open the approprA_countte data file */
   input_fs = openDataFile();   
   
   
   strain_offset = get_snp_offset();
   
   snp.iEncodedSNP = malloc(sizeof(int) * g_data_info.num_strains);
   
   /* open output file */
   if (prefix != NULL)
   {
      strncpy(path, prefix, MAXPATHLEN);
      strncpy(path_filtered, prefix, MAXPATHLEN);
   }
   
   /* name the file based on which path we used to fill the data */
   if (path_option == PATH_VITERBI)
   {
      snprintf(path+strlen(path), MAXPATHLEN,"viterbi_filled");
      snprintf(path_filtered+strlen(path_filtered), MAXPATHLEN,"viterbi_filled");
   }
   else
   {
      snprintf(path+strlen(path), MAXPATHLEN,"ms_filled");
      snprintf(path_filtered+strlen(path_filtered), MAXPATHLEN,"ms_filled");
   }
   
   if (sort_option == NO_SORT)
   {
      snprintf(path+strlen(path), MAXPATHLEN,".csv");
      snprintf(path_filtered+strlen(path_filtered), MAXPATHLEN,"_filtered.csv");
   }
   else
   {
      snprintf(path+strlen(path), MAXPATHLEN,"_sorted.csv");
      snprintf(path_filtered+strlen(path_filtered), MAXPATHLEN,"_sorted_filtered.csv");
   }
   
   printf("writing filled Data to: \n   %s\n   %s\n", path, path_filtered);
   
   filled_output_fs = fopen(path, "w");
   
   if (filled_output_fs == NULL)
   {
      fprintf(stderr, "Unable to open filled output file:\n\t%s\n\t%s\n", 
              path, strerror(errno));
      exit(errno);
   }
   
   filtered_output_fs = fopen(path_filtered, "w");
   
   if (filtered_output_fs == NULL)
   {
      fprintf(stderr, "Unable to open filled output file:\n\t%s\n\t%s\n", 
              path_filtered, strerror(errno));
      exit(errno);
   }
   
   if (prefix != NULL)
   {
      strncpy(path, prefix, MAXPATHLEN);
   }
   if (path_option == PATH_VITERBI)
   {
      snprintf(path+strlen(path), MAXPATHLEN,"viterbi_filling_prob1");
   }
   else
   {
      snprintf(path+strlen(path), MAXPATHLEN,"ms_filling_prob1");
   }  
   
   if (sort_option == NO_SORT)
   {
      snprintf(path+strlen(path), MAXPATHLEN,".csv");
   }
   else
   {
      snprintf(path+strlen(path), MAXPATHLEN,"_sorted.csv");
   }

   filling_prob1_fs = fopen(path, "w");
   
   if (filling_prob1_fs == NULL)
   {
      fprintf(stderr, "Unable to open filling probability output file:\n\t%s\n\t%s\n", 
              path, strerror(errno));
      exit(errno);
   }
   

   if (prefix != NULL)
   {
      strncpy(path, prefix, MAXPATHLEN);
   }
   
   if (path_option == PATH_VITERBI)
   {
      snprintf(path+strlen(path), MAXPATHLEN,"viterbi_filling_prob2");
   }
   else
   {
      snprintf(path+strlen(path), MAXPATHLEN,"ms_filling_prob2");
   } 
   
   if (sort_option == NO_SORT)
   {
      snprintf(path+strlen(path), MAXPATHLEN,".csv");
   }
   else
   {
      snprintf(path+strlen(path), MAXPATHLEN,"_sorted.csv");
   }

   filling_prob2_fs = fopen(path, "w");
   if (filling_prob2_fs == NULL)
   {
      fprintf(stderr, "Unable to open filling probability output file:\n\t%s\n\t%s\n", 
              path, strerror(errno));
      exit(errno);
   }
   
   
   

   /* copy the header */

   line_buffer = csvGetLine(input_fs);
   snp_id_column = csvGetID_Column();
   for (int i = 0; i < csvNumFields(); i++)
   {
      fprintf(filled_output_fs, "%s", csvField(i));
      fprintf(filtered_output_fs, "%s", csvField(i));
      fprintf(filling_prob1_fs, "%s", csvField(i));
      fprintf(filling_prob2_fs, "%s", csvField(i));

      if (i < csvNumFields() - 1)
      {
         fprintf(filled_output_fs, ",");
         fprintf(filtered_output_fs, ",");
         fprintf(filling_prob1_fs, ",");
         fprintf(filling_prob2_fs, ",");

      }
      else
      {
         fprintf(filled_output_fs, "\n");
         fprintf(filtered_output_fs, "\n");
         fprintf(filling_prob1_fs, "\n");
         fprintf(filling_prob2_fs, "\n");  
      }
   }

   chr_column = csvGetChrColumn();
   
   
   /* read data, fill, write out */

 
   /* read the data line by line */
   index = 0;
   while((line_buffer = csvGetLine(input_fs)) != NULL && good_snp_idx < good_snp_count)
   {  
      
      /* if the chromosome is not right skip */
      if (getChrNum(csvField(chr_column)) != chr)
      {
         index++;
         continue;
      }
      
      /* tokenize string, using ',' as a delimiter */
      
      out_line[0] = '\0';
      strain_index = 0;
      for (int i = 0; i < csvNumFields(); i++)
      {
         
         /* copy non-sequence columns */
         if (i < strain_offset - 1)
         {
            strcat(out_line, csvField(i));
            strcat(out_line, ",");
         }
         else 
         {
            char *word = csvField(i);
            /* else copy the unfilled SNP */
            snp_str[strain_index++] = word[0];   
         }
      }
       
      snp_str[strain_index] = '\0';
      
      fprintf(filling_prob1_fs, "%s", out_line);
      fprintf(filling_prob2_fs, "%s", out_line);


      if (snps[good_snp_idx].index == index)
      {
         psnp = &snps[good_snp_idx];

         /* fill in missing information in the SNP */   
         fillSNP(snp_str, &snps[good_snp_idx], g_data_info.num_strains);

         for (int iStrain = 0; iStrain < g_data_info.num_strains; iStrain++)
         {
            fprintf(filling_prob1_fs, "%f", fill_probabilities[good_snp_idx][iStrain]);

            if (islower(snp_str[iStrain]))
            {
               fprintf(filling_prob2_fs, "%f", 
                       fill_probabilities[good_snp_idx][iStrain]);
            }
            else
            {
               fprintf(filling_prob2_fs, "-");
            }
               
            if (iStrain == g_data_info.num_strains - 1)
            {
               fprintf(filling_prob1_fs, "\n");
               fprintf(filling_prob2_fs, "\n");
            }
            else
            {
               fprintf(filling_prob1_fs, ",");
               fprintf(filling_prob2_fs, ",");
            }
               
         }
            
      }
      else
      {
         /* snp was removed, we need to read it in so we can fill it... */
                              
         encodeSNP(snp_str, &snp, g_data_info.filter_flags);
         psnp = &snp;
         for (int iStrain = 0; iStrain < g_data_info.num_strains; iStrain++)
         {
            fprintf(filling_prob1_fs, "-");
            fprintf(filling_prob2_fs, "-");            
     
            if (iStrain == g_data_info.num_strains - 1)
            {
               fprintf(filling_prob1_fs, "\n");
               fprintf(filling_prob2_fs, "\n");
            }
            else
            {
               fprintf(filling_prob1_fs, ",");
               fprintf(filling_prob2_fs, ",");
            } 
            
         }

      }
         
      /* this is needed for RIL cases where one of the parents was N, but
         we could determine the Allele from the other strains. Since 
         that parent is now encoded as 1 or 0, it won't get filled 
         in fillSNP, but we need to fill it with the determined 
         Allele, not keep the N */
      if (g_data_info.num_parent_strains == 2 && 
          psnp->iStatus == SNP_GOOD)
      {
         if (toupper(snp_str[0]) == 'N' && psnp->cAlleleArray[0] != 'N')
         {
            snp_str[0] = tolower(psnp->cAlleleArray[0]);
         }
         if (toupper(snp_str[1]) == 'N'&& psnp->cAlleleArray[1] != 'N')
         {
            snp_str[1] = tolower(psnp->cAlleleArray[1]);
         }
      }

         
      /*  copy SNP into our output buffer, adding in the commas */
      strcpy(out_line_filtered, out_line);
      for(int j = 0; j < g_data_info.num_strains; j++)
      {
         sprintf(&out_line[strlen(out_line)], "%c,", snp_str[j]);
         if (psnp->iStatus != SNP_GOOD || 
             isupper(snp_str[j]) ||
             (psnp->iStatus == SNP_GOOD && 
             fill_probabilities[good_snp_idx][j] >= threshold) ) {
            sprintf(&out_line_filtered[strlen(out_line_filtered)], "%c,", snp_str[j]);
         }
         else
         {
            sprintf(&out_line_filtered[strlen(out_line_filtered)], "N,");
         }
           
      }
      out_line[strlen(out_line) - 1] = '\n';
      out_line_filtered[strlen(out_line_filtered) - 1] = '\n';

      fputs(out_line, filled_output_fs);  
      fputs(out_line_filtered, filtered_output_fs);  

         
      if (psnp->iStatus == SNP_GOOD)
      {
         good_snp_idx++;    
      }

      index++;
       
   }

   fclose(input_fs);
   fclose(filtered_output_fs);
   fclose(filled_output_fs);
   fclose(filling_prob1_fs);
   fclose(filling_prob2_fs);
   
}


/*******************************************************************************
 * writePath - write the haplotype path out to a file
 *
 * INPUTS: snpt_t snps[] - SNP data
 *         int **path_matrix - 2D array containing the path for each sequence
 *         int num_snps - number of SNPs
 *         int num_strains - number of strains
 *         const char prefix[] - prefix for output file name
 *         int path_option - was the path computed with viterbi or maxsmth
 *         int chr - which chromosome this pathis for
 *         int sort_option - sorting option, used to generate file name
 *         
 *
 * RETURNS: void
 *
 * ASSUMES: that we can open up the original data file to reconstruct the 
 *    annotation columns
 *
 * EFFECTS: writes path info to file. may overwrite older path file
 *
 * ERROR CONDITIONS: will exit if unable to open file
 *
 * COMMENTS:
 *
 ******************************************************************************/
void writePath(snp_t *snps, int **path_matrix, int num_snps, 
               int num_strains, const char prefix[], int path_option, 
               int chr, int sort_option, double ***emission_matrix)
{
   char path[MAXPATHLEN];
   char *line_buffer;
   char out_line[LINE_BUFFER_SIZE];
   
   FILE *filled_output_fs; 
   FILE *input_fs;
   
   int strain_offset;
   int chr_column = 0;
   int haplotypes;
   int emission_types;
   int snp_index = 0;
   int path_index = 0;
   char *fill_rule;
   int index;
   
   /* open the data file so we can get the annotation columns*/

   input_fs = openDataFile();   
   strain_offset = get_snp_offset();

   emission_types = get_num_emission_types();
   haplotypes = get_num_haplotypes();
   fill_rule = malloc((haplotypes + 1) * sizeof(char));
   
   /* open up the output file */
   if (prefix != NULL)
   {
      strncpy(path, prefix, MAXPATHLEN);
   }
   
   if (path_option == PATH_VITERBI)
   {
      snprintf(path+strlen(path), MAXPATHLEN,"viterbi_path");
   }
   else
   {
      snprintf(path+strlen(path), MAXPATHLEN,"ms_path");
   }
   
   if (sort_option == NO_SORT)
   {
      snprintf(path+strlen(path), MAXPATHLEN,".csv");
   }
   else
   {
      snprintf(path+strlen(path), MAXPATHLEN,"_sorted.csv");
   }
   
   printf("writing path to: %s\n", path);
      
   if ( (filled_output_fs = fopen(path, "w")) == NULL)
   {
      /* TODO we probably shouldn't just bail out like this... */
      printf("Unable to open path file:\n\t%s\n\t%s\n", path, strerror(errno));
      exit(errno);
   }
   
   /* copy the header*/

   line_buffer = csvGetLine(input_fs);
   chr_column = csvGetChrColumn();
      
   out_line[0] = 0;    
      
   for (int i = 0; i < strain_offset - 1; i++)
   {
      fprintf(filled_output_fs, "%s,", csvField(i));
   }
   fprintf(filled_output_fs, "#Haplotypes,filling.rule,");
   for (int i = strain_offset - 1; i < csvNumFields(); i++)
   {
      fprintf(filled_output_fs, "%s", csvField(i));
      if (i < csvNumFields() - 1)
      {
         fprintf(filled_output_fs, ",");
      }
      else
      {
         fprintf(filled_output_fs, "\n");
      }
   }

   
   /* read data, write out path with annotation data*/
   index = 0;
   while((line_buffer = csvGetLine(input_fs)) != NULL)
   {
      /* if not the right chromosome then we skip
         the SNP */
      if (getChrNum(csvField(chr_column)) != chr)
      {
         index++;
         continue;
      }


      out_line[0] = 0;
      for (int i = 0; i < strain_offset - 1; i++)
      {
         sprintf(out_line+strlen(out_line), "%s,", csvField(i));   
      }
                  
      if (snps[snp_index].index == index)
      {
         /* modifying code to no longer store "bad" SNPs in the array... */      
         assert(snps[snp_index].iStatus == SNP_GOOD);
        
         calculate_fill_rule(snps[snp_index], emission_matrix[snp_index+1], 
                             haplotypes, emission_types, fill_rule);
        
         /* copy number of haplotypes and fillin.rule into buffer */
         sprintf(&out_line[strlen(out_line)], "%d,%s", 
                 snps[snp_index].iNumHaplotypes, 
                 fill_rule);
         
         /*  copy inferred haplotypes into output buffer,
             adding in the commas */
         for(int j = 0; j < g_data_info.num_strains; j++)
         {

            sprintf(&out_line[strlen(out_line)], ",%d", 
                    path_matrix[path_index][j]);
         }
         path_index++;
         snp_index++;
      }
      else
      {
         /* this snp was removed, so we are just puting in a placeholder.
            use '-' for the # haplotypes and all strain inferred haplotypes
            (g_data_info.num_strains + 1 total columns) */
         strcat(out_line, "-");
         for (int j = 0; j < g_data_info.num_strains + 1; j++)
         {
            strcat(out_line, ",-");
         }
      }
      
      /* add a newline */
      out_line[strlen(out_line)  + 1] = '\0';
      out_line[strlen(out_line)] = '\n';

      fputs(out_line, filled_output_fs);  
      index++;
      }
  
   free(fill_rule);
   fclose (input_fs);
   fclose (filled_output_fs);
}

/*******************************************************************************
 * writeLambdaOutpBinary - write the lambda/outp values into a binary file
 *
 * INPUTS:  char path[] - path to write data to
 *          double ***transition_matrix - lambda values
 *          double ***emission_matrix - outp values
 *          int num_haplotypes - number of haplotypes
 *          int num_snps - number of SNPs
 *          int iEmissionTypes = number of emission types (2 for random missing,
 *             3 for missing as emission)
 *
 * RETURNS: void
 *
 * ASSUMES: transition_matrix, emission_matrix are allocated
 *
 * EFFECTS: creates or overwrites file to save binary data
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void writeLambdaOutpBinary(char path[], double ***transition_matrix, 
                           double ***emission_matrix,  
                           int num_haplotypes, int num_snps, 
                           int emission_types)
{
   FILE* output_fs;
   
   output_fs = fopen(path, "wb");
   
   
   if (output_fs == NULL)
   {
      perror("Unable to open Lambda output file");
      exit(errno);
   }
   
   
   /* beginning */
   
   fwrite(&transition_matrix[0][0][1], sizeof(double), num_haplotypes, output_fs);
   
   fwrite(&emission_matrix[0][0][0], sizeof(double), emission_types, output_fs);
         
   
   
   
   /* do SNP1 - SNPT-1 */
   for (int i = 1; i < num_snps; i++)
   {
      for (int j = 0; j < num_haplotypes; j++)
      {
         
         fwrite(&transition_matrix[i][j][0], sizeof(double), num_haplotypes, 
                output_fs);
         
         fwrite(&emission_matrix[i][j][0], sizeof(double), emission_types, 
                output_fs);

      }
   }
   /* do SNP T */
   
   for (int i = 0; i < num_haplotypes; i++)
   {
      fwrite(&transition_matrix[num_snps][i][0], sizeof(double), 1, output_fs);

      fwrite(&emission_matrix[num_snps][i][0], sizeof(double), emission_types, 
             output_fs);
                  
   }
   
   /* do end */
   
   fwrite(&emission_matrix[num_snps+1][0][0], sizeof(double), emission_types, 
          output_fs);
         
   
   fclose (output_fs);


}

/*******************************************************************************
 * readLambdaOutpBinary - read in a binary file with lambda / outp values
 *
 * INPUTS:  char path[] - path to binary file
 *          double ***transition_matrix
 *          double ***emission_matrix
 *          int num_haplotypes = number of haplotypes to read in
 *          int num_snps - number of SNPs to read in
 *          int emission_types - number of emission types in outp
 *
 * RETURNS: void
 *
 * ASSUMES:
 *
 * EFFECTS: modifies contents of transition_matrix and emission_matrix
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void readLambdaOutpBinary(char path[], double ***transition_matrix, 
                          double ***emission_matrix,  
                          int num_haplotypes, int num_snps, 
                          int emission_types)
{
   FILE* input_fs;
   
   input_fs = fopen(path, "rb");
   
   if (input_fs == NULL)
   {
      perror("Unable to open Lambda input file");
      exit(errno);
   }
   
   
   /* beginning */
   
   fread(&transition_matrix[0][0][1], sizeof(double), num_haplotypes, input_fs);
   
   fread(&emission_matrix[0][0][0], sizeof(double), emission_types, input_fs);
         
   
   
   
   /* do SNP1 - SNPT-1 */
   for (int i = 1; i < num_snps; i++)
   {
      for (int j = 0; j < num_haplotypes; j++)
      {
         
         fread(&transition_matrix[i][j][0], sizeof(double), num_haplotypes, 
               input_fs);
         
         fread(&emission_matrix[i][j][0], sizeof(double), emission_types, 
               input_fs);

      }
   }
   /* do SNP T */
   
   for (int i = 0; i < num_haplotypes; i++)
   {
      fread(&transition_matrix[num_snps][i][0], sizeof(double), 1, input_fs);

      fread(&emission_matrix[num_snps][i][0], sizeof(double), emission_types, 
            input_fs);
                  
   }
   
   /* do end */
   
   fread(&emission_matrix[num_snps+1][0][0], sizeof(double), emission_types, 
         input_fs);
         
   
   fclose (input_fs);


}

/**** I N T E R N A L    F U N C T I O N S *****/



static int check_snp(char *snp_str, int filter_flags)
{
   char allele_A = '\0';
   char allele_B = '\0';
   int A_count = 0;
   int B_count = 0;
   int num_missing  = 0;
   double missingness;
   char *tmp;
   
   /* this function is desctructive, so we make a copy of the string first */
   tmp = strdup(snp_str);
   
   for (int i = 0; i < g_data_info.num_strains; i++)
   {
      /* turn anything that isn't AGCT to N */
      if (toupper(tmp[i]) != 'A' && toupper(tmp[i]) != 'G' 
          &&  toupper(tmp[i]) != 'C' && toupper(tmp[i]) != 'T')
      {
         tmp[i] = 'N';
      }

     
      if (toupper(tmp[i]) == 'N')
      {
         num_missing++;
      }
      else if (isTraining(i) && A_count == 0)
      {
         allele_A = toupper(tmp[i]);
         A_count++;
      }
      else if (isTraining(i) && tmp[i] == toupper(allele_A))
      {
         A_count++;
      }
      else if (isTraining(i) && tmp[i] != allele_A && B_count == 0)
      {
         allele_B = toupper(tmp[i]);
         B_count++;
      }
      else if (isTraining(i) && tmp[i] == toupper(allele_B))
      {
         B_count++;
      }

      
      
   }
   
   
   /* make sure there are not more than two non-N alleles */
   

   if (A_count + B_count + num_missing < g_data_info.num_strains)
   {
      free(tmp);
      return SNP_TOO_MANY_NT;
   }

   
   if (A_count == 0 && B_count == 0 && g_data_info.num_parent_strains == 0)
   {
      free(tmp);
      return SNP_TOO_MANY_MISSING;
   }
   
 
   /* missingness = num missing / num strains */
   missingness =  (double)num_missing / g_data_info.num_strains;
   
   
   if (missingness > g_data_info.max_missingness)
   {
      free(tmp);
      return SNP_TOO_MANY_MISSING;
   }
   
   
   if (g_data_info.num_parent_strains == 0)
   {
      
      /* only one non-N nucleotide, non RIL, and we are removing constant SNPs */
      
      if (B_count == 0 && ((filter_flags & KEEP_CONSTANT_FLAG) == 0))
      {
         free(tmp);
         return SNP_ONE_NT; 
      }

   }
   else  /* RIL SNPs */
   {  
   
      /* by convention the first two strains are the parents */
      char c0 = toupper(tmp[PARENT_1]);
      char c1 = toupper(tmp[PARENT_2]);
   
      if (c0 == 'N' && c1 == 'N')
      {
         free(tmp);
         return SNP_RIL_MISSING_PARENTS;
      }
      else if ( c0 == c1)
      {
         free(tmp);
         return SNP_RIL_CONSTANT;  
      }
      else if (c0 == 'N' || c1 == 'N') {
         free(tmp);
         return SNP_RIL_MISSING_PARENT;
      }
      
      for (int i = 0; i < g_data_info.num_strains; i++)
      {
         /* make sure all non-N values match one of the parents */
         if (toupper(tmp[i]) != 'N' && toupper(tmp[i]) != c0 && 
             toupper(tmp[i]) != c1)
         {
            free(tmp);
            return SNP_RIL_NOMATCH;
         }
      }
   }
   
   /* if SNP is good then this must be true or we have a logic error */
   assert(A_count + B_count + num_missing == g_data_info.num_strains);
   
   free(tmp);
   return SNP_GOOD;
}


/*******************************************************************************
 * encodeSNP - takes a SNP and encodes it in the integer format
 *   0 - major allele, 1 - minor allele, 2 - missing
 *
 * INPUTS: char snp_str[] - SNP string
 *         snpt_t snpEnc - SNP encoded into struct
 *         int iFiltlerFlags - various flags that alter filtering behavior
 *         cjar *snp_id - the SNP id if available, otherwise "NA"
 *
 * RETURNS: void
 *
 * ASSUMES: snp_str contains a valid SNP
 *
 * EFFECTS: alters psnpEnc
 *
 * ERROR CONDITIONS: if SNP contains more than two nucleotide then the function 
 *   will print an error message and cause the program to exit
 *
 * COMMENTS:
 *
 ******************************************************************************/
static void encodeSNP(char snp_str[], snp_t *snp_enc, int filter_flags)
{
   /* first we need to find the '0' and '1' values */
   
   char allele_A = '\0';
   char allele_B = '\0';
   int A_count = 0;
   int B_count = 0;
   int num_missing  = 0;
   double missingness;
   char c0, c1;
   
   snp_enc->iStatus = check_snp(snp_str, filter_flags);
   

 

   /* we know if the SNP is "good" or not, but we aren't quite ready to return
      even if the SNP is bad...*/

   
   /* scan over the SNP, setting all missing to 2, and count the two alleles */
   for (int i = 0; i < g_data_info.num_strains; i++)
   {
      /* turn anything that isn't AGCT to N */
      if (toupper(snp_str[i]) != 'A' && toupper(snp_str[i]) != 'G' 
          &&  toupper(snp_str[i]) != 'C' && toupper(snp_str[i]) != 'T')
      {
         snp_str[i] = 'N';
      }
          
      if (toupper(snp_str[i]) == 'N')
      {
         snp_enc->iEncodedSNP[i] = 2;
         num_missing++;
      }
      else if (isTraining(i) && A_count == 0)
      {
         allele_A = toupper(snp_str[i]);
         A_count++;
      }
      else if (isTraining(i) && snp_str[i] == toupper(allele_A))
      {
         A_count++;
      }
      else if (isTraining(i) && snp_str[i] != allele_A && B_count == 0)
      {
         allele_B = toupper(snp_str[i]);
         B_count++;
      }
      else if (isTraining(i) && snp_str[i] == toupper(allele_B))
      {
         B_count++;
      }
   
      
   }
   

   /*
      if the status is not SNP_GOOD we will return, but for a few different 
      statuses we may need to set some values used by the filling process first
    */
   if (snp_enc->iStatus != SNP_GOOD)
   {
   
      /* this is so our filling function will work for a constant SNP */
      if (snp_enc->iStatus == SNP_ONE_NT)
      {
         snp_enc->cAlleleArray[0] = allele_A;
      }
      else if (snp_enc->iStatus == SNP_RIL_MISSING_PARENTS)
      {
         snp_enc->cAlleleArray[0] = 'N';
         snp_enc->cAlleleArray[1] = 'N';
      }
      else if (snp_enc->iStatus == SNP_RIL_MISSING_PARENT ||
               snp_enc->iStatus == SNP_RIL_NOMATCH)
      {
         snp_enc->cAlleleArray[0] = snp_str[PARENT_1];
         snp_enc->cAlleleArray[1] = snp_str[PARENT_2];
      }
      else if (snp_enc->iStatus == SNP_RIL_CONSTANT) {
         /* both alleles are the same, we only use cAlleleArray[0] when filling */
         snp_enc->cAlleleArray[0] = toupper(snp_str[PARENT_1]);
      }
      
      /* no need to encode the SNP, we won't be using it in the model */
      return;
   }
  
 
   /* missingness = num missing / num strains */
   missingness =  (double)num_missing / g_data_info.num_strains;

   
   
   /* now we can find which is the major and which is the minor allele. */

   
   if (g_data_info.num_parent_strains == 0)
   {
      
      /* only one non-N nucleotide, non RIL */
      
      if (B_count == 0)
      {   

         /* this is a "constant" SNP, but if we end up with any 
            "1"s in the encoded SNP we want to fill them in with 'n' since we 
            don't know if the SNP is really constant or we are just missing 
            all information for strains with the second allele */
         allele_B = 'N';

      }

      
      /* if SNP is constant then A_count is always bigger than B_count, so the major
         allele will be the one observed allele, the minor allele will be 'N'
         ...  If the couts A_count and B_count are the same then we order major/minor 
         using alphabetical ordering */
      if (A_count > B_count || (A_count == B_count && allele_A < allele_B))
      {
         c0 = allele_A;
         c1 = allele_B;
      }      
      else
      {
         c0 = allele_B;
         c1 = allele_A;
      }
      
      /* the allele array will be used in the filling process later */
      snp_enc->cAlleleArray[0] = c0;
      snp_enc->cAlleleArray[1] = c1;

      
   }
   else 
   {  
      /* this is RIL data */
   
      c0 = toupper(snp_str[PARENT_1]);
      c1 = toupper(snp_str[PARENT_2]);
   
      snp_enc->cAlleleArray[0] = c0;
      snp_enc->cAlleleArray[1] = c1;
   
   
      /* TODO make "guessing" missing parent genotype a user specified option */
      if (c0 == 'N')
      {
         /* since we remove SNPs with poth parents missing, this must be true */
         assert(c1 != 'N');
      
         /* if parent 1 type is missing, try to guess by looking at other
            genotypes */
         for (int i = PARENT_2 + 1; i < g_data_info.num_strains; i++)
         {
            if (toupper(snp_str[i]) != c1 && toupper(snp_str[i]) != 'N')
            {
               c0 = toupper(snp_str[i]);
               snp_enc->cAlleleArray[0] = c0;
               snp_str[PARENT_1] = tolower(c0);
               break;
            }
         }
        
      }
      else if (c1 == 'N')
      {
         /* since we remove SNPs with poth parents missing, this must be true */
         assert(c0 != 'N');
        
         /* if parent type 2 is missing try to guess by looking at other 
           genotypes */
         for (int i = PARENT_1 + 1; i < g_data_info.num_strains; i++)
         {
            if (toupper(snp_str[i]) != c0 && toupper(snp_str[i]) != 'N')
            {
               c1 = toupper(snp_str[i]);
               snp_enc->cAlleleArray[1] = c1;
               snp_str[PARENT_2] = tolower(c1);
               break;
            }
         }       
      }
      
   }
   
   

   /* encode major allele as 0, minor allele as 1. leave N as 2 */
   for (int i = 0; i < g_data_info.num_strains; i++)
   {
      /* if c0 == c1, then we encode as 0 */
      if (toupper(snp_str[i]) == toupper(c0) && c0 != 'N')
      {
         snp_enc->iEncodedSNP[i] = 0;
      }
      else if (toupper(snp_str[i]) == toupper(c1) && c1 != 'N')
      {
         snp_enc->iEncodedSNP[i] = 1;
      }
      else if (snp_str[i] == 'N' || snp_str[i] == 'n')
      {
         snp_enc->iEncodedSNP[i] = 2;
      } 
     
   } 

}



/*******************************************************************************
 * initNumStrains - returns number of strains in a single .csv datafile
 *
 * INPUTS: int strain_offset - column offset of first sequence column
 *
 * RETURNS: integer value of number of strains in the data set
 *
 * ASSUMES: ddataset has been initA_countlized and is valid
 *
 * EFFECTS:none
 *
 * ERROR CONDITIONS: will exit if unable to open data file
 *
 * COMMENTS: not to be called directly - this is called by getNumStrains
 *
 ******************************************************************************/
static int initNumStrains(int strain_offset)
{
     
   char *str;
   g_data_info.num_strains = 0;
   FILE* input_fs = NULL;
      
   if ( (input_fs = fopen(g_data_info.data_path, "r")) == NULL)
   {
      printf("Unable to open data file:\n\t%s\n\t%s\n", g_data_info.data_path,
             strerror(errno));
      exit(errno);
   }
         
   if ((str = csvGetLine(input_fs)) != NULL)
   {
      g_data_info.num_strains = csvNumFields() - (strain_offset - 1);
   }
   else 
   {
      fprintf(stderr, "Error calculating number of strains in input file\n");
      exit(EHMM);
   }
   
   return g_data_info.num_strains;
}



/*******************************************************************************
 * initStrainNames - get the strain names from a single .csv data file
 *
 * INPUTS: int strain_offset - colunm offset for first sequence
 *
 * RETURNS: void
 *
 * ASSUMES: first line in the file is a properly formatted header
 *
 * EFFECTS: allocates memory for each strain name string, copies pointers 
 *    to these strings into the g_data_info struct
 *
 * ERROR CONDITIONS: will exit if unable to open the datafile
 *
 * COMMENTS: this is only called in initA_countlizeDataSource
 *
 ******************************************************************************/
static void initStrainNames(int strain_offset)
{

   char *line_buffer;
   int strain = 0;
   /* input data is a single file. Get the strain names*/
   FILE* input_fs = NULL;
   
   if ( (input_fs = fopen(g_data_info.data_path, "r")) == NULL)
   {
      printf("Unable to open data file:\n\t%s\n\t%s\n", g_data_info.data_path,
             strerror(errno));
      exit(errno);
   }
      
   /* read the first line out of the file */   
   if ((line_buffer = csvGetLine(input_fs)) != NULL )
   {
      for (int i = 0; i < csvNumFields(); i++)
      {
         if (i >= strain_offset - 1)
         {
            g_data_info.strain_names[strain++] = strdup(csvField(i));
         }
      }

   }
   fclose (input_fs);
}



/*******************************************************************************
 * getChromosomeSizes - get the number of SNPs in each chromosome for a 
 *   single .csv data file
 *
 * INPUTS: void
 *
 * RETURNS: void
 *
 * ASSUMES: data file is a valid format
 *
 * EFFECTS: copies sizes of each chromosome into g_data_info struct
 *
 * ERROR CONDITIONS: will cause program to exit if unable to open file
 *
 * COMMENTS: this is only called in initA_countlizeDataSource
 *
 ******************************************************************************/
static void initChromosomeSizes()
{
   int chr_column = 0;
   char *line_buffer;
   /* input file is a unified file. */
   FILE* input_fs = NULL;

   /* get sizes of chromosomes for single file data format */
   
   
   for (int i = 0; i < NUM_CHROMOSOMES; i++)
      g_data_info.chr_sizes[i] = 0;

      
   if ( (input_fs = fopen(g_data_info.data_path, "r")) == NULL)
   {
      printf("Unable to open data file:\n\t%s\n\t%s\n", g_data_info.data_path,
             strerror(errno));
      exit(errno);
   }
         
   int line_num = 0;
   while ((line_buffer = csvGetLine(input_fs)) != NULL )
   {  
         
      if (line_num++ == 0)
      {
         chr_column = csvGetChrColumn();
       
      } 
      else
      {
         int index = getChrNum(csvField(chr_column));
         if (index >= 0)
         {
            g_data_info.chr_sizes[index-1]++;
         }
         else 
         {
            printf("Error parsing chromosome field %s\n", csvField(chr_column));
            exit(EHMM);
         }
         
      }
   }

}


/*******************************************************************************
 * openDataFile - open single .csv data file
 *
 * INPUTS: none  
 *
 * RETURNS: open file stream (for reading)
 *
 * ASSUMES: datasource has been validated 
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS: will cause program to exit if file can not be opened.
 *
 * COMMENTS: simple function, but provided to complement openDataDir
 *
 ******************************************************************************/
static FILE *openDataFile()
{
   FILE *input_fs = NULL;
   
   if ( (input_fs = fopen(g_data_info.data_path, "r")) == NULL)
   {
      printf("Unable to open data file:\n\t%s\n\t%s\n", g_data_info.data_path,
             strerror(errno));
      exit(errno);
   }
   
   return input_fs;
}
