/*
   File Name: hmmUtil.c
   System: hmmSNP
   Programmer: Glen Beane with assistance from Jin Szatkeiwicz
   Date developed:
   Purpose:
   Overview: this file contains various helper functions used throughout
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


#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>

#include "constants.h"
#include "hmmUtil.h"

/* string utility functions */


/*******************************************************************************
 * int isDir(const char *path)
 *
 * INPUTS:  const char *path: path to file we want to test
 *
 * RETURNS: TRUE if path points to a directory, FALSE otherwise
 *          ELSTAT if call to stat() fails
 *
 * ASSUMES: path is valid input for stat()
 *
 * EFFECTS: may write error message to stderr
 *
 * ERROR CONDITIONS: will return error code if stat() function fails
 *
 * COMMENTS: TRUE and FALSE are defined in constants.h
 *
 ******************************************************************************/
int isDir(const char *path)
{
   struct stat sbuf;

   if (path == NULL) 
   {
      return FALSE;
   }
  
   if (stat(path, &sbuf) == -1) 
   {
      fprintf(stderr, "stat() Failed.\n");
      return ELSTAT;
   }

   if(S_ISDIR(sbuf.st_mode))
   {
      return TRUE;
   }
   return FALSE;
} 




/*******************************************************************************
 * int getChrNum(const char *chr_str)
 *
 * INPUTS: char chr_str[] - chromosome field from a line in the input file
 *         
 *
 * RETURNS: integer value of chromosome number
 *
 * ASSUMES:
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS: returns a chromosome number from a csvField
 ******************************************************************************/
int getChrNum(char *chr_str)
{

   int chr = 0;
   int i;
   
   char *tmp = strdup(chr_str); /* make a copy - we'll be modifying it */
   char *p;
 
   /* get pointer to tmp because we'll be moving the pointer to skip extra 
      stuff at the beginning of the field, we'll need the original to 
      free it later */
   p = tmp;            

   i = strlen(p) - 1;
   while (!isdigit(p[i]) && !isalpha(p[i]))
   {
      p[i--] = '\0';
      if (i < 0)
      {
         return -1;
      }
   }

   /* skip extra LEADING characters.  chromosome can contain a number of 
      extra characters (usually Chr) followed by the chromosome name. 
      The chromosome name must occur last in the string. */
   
   
   while (!isdigit(p[0]) && strlen(p) > 1)
   {
      p++;
   }
      
   if ( !isdigit(p[0]) )
   {
      if (tolower(p[0]) == 'x')
      {
         chr = CHR_X;
      }
      else if (tolower(p[0]) == 'y')
      {
         chr = CHR_Y;
      }
      else if (tolower(p[0]) == 'm')
      {
         chr = CHR_M;
      }
      else
      {
         /* bad chr string */
         chr = -1;
      }
   }
   else
   {
      char *invalid_char;
      chr = (int)strtol(p, &invalid_char, 10);
      if (*invalid_char != '\0')
      {  
         chr = -1;
      }
   }
         
   
   free(tmp);
   return chr;

}
/*******************************************************************************
 * getSNP(const char *line, char snp_str[], int iOffsset, int strain_flags[])
 *
 * INPUTS: const char *line - string containin a line from the data file
 *         char snp_str - array to copy SNP into
 *         int offset - offset of first sequence in data set
 *         int strain_flags[] - array tells if we are using each sequence
 *
 * RETURNS: void.  returns SNP via snp_str[]
 *
 * ASSUMES:
 *
 * EFFECTS: alters contents snp_str[]
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS: this function takes a complete line from a datafile and pulls out 
 *  the SNP data
 *
 ******************************************************************************/
void getSNP(const char *line, char snp_str[], int offset, 
            int strain_flags[])
{
   char *tmp_str = strdup(line);
   int index = 0;
   char *token, *idx;
   idx = tmp_str;

      
   token = strtok(idx, ",");
 
   for (int i = 1; i < offset; i++)
   {
      token = strtok(idx,",");
   }

   /* grab SNP data */
   int i = 0;
   while (token != NULL)
   {
      if (strain_flags[i++] == TRUE)
      {
         snp_str[index++] = token[0];
      }
      token = strtok(idx,",");
   }
         

}


/*******************************************************************************
 * setupOutputPrefix - sets up a string with the approprate output file prefix
 *
 * INPUTS:   int chr - chromosome number
 *           char user_prefix - user supplied prefix (may be NULL)
 *
 * RETURNS:  pointer to prefix string
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
char *setupOutputPrefix(int chr, char user_prefix[])
{
   char *prefix; 
   
   if (user_prefix != NULL)
   {
      prefix = user_prefix;
   }
   else
   {
      prefix = malloc(sizeof(char) * (MAX_PREFIX + 1));
      
      if (prefix == NULL)
      {
         fprintf(stderr, "Unable to allocate Prefix\n");
         exit (EHMM_NOMEM);
      }
      switch (chr)
      {
         case CHR_X:
            strcpy(prefix, "chrX_");         
            break;
         case CHR_Y:
            strcpy(prefix, "chrY_");
            break;
         case CHR_M:
            strcpy(prefix, "M_");
            break;
         default:
            snprintf(prefix, MAX_PREFIX, "chr%d_", chr);
      }   
   }
   
   return prefix;
}



/* memory utilities */


/*******************************************************************************
 * allocate2Dd - allocate 2D array of doubles
 *
 * INPUTS: int height - number of rows
 *         int width - number of columns
 *
 * RETURNS: pointer to allocated memory
 *
 * ASSUMES:
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS: returns NULL if unable to allocate memory.
 *
 * COMMENTS:
 *
 ******************************************************************************/   
double **allocate2Dd(int height, int width)
{

   double **p = malloc(sizeof(double*) * height);
   
   if (p == NULL)
      return NULL;
   
   for (int i = 0; i < height; i++)
   {
      p[i] = malloc(sizeof(double) * width);
      if (p[i] == NULL)
      {
         for (int j = 0; j < i; j++)
         {
            free(p[j]);
         }
         return NULL;
      }
   }
   
   return p;

}

/*******************************************************************************
 * free2Dd - free a 2D double array
 *
 * INPUTS: double **p - matrix to free
 *         int height - number of rows
 *
 * RETURNS: void
 *
 * ASSUMES:
 *
 * EFFECTS:  frees memory
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void free2Dd(double **p, int height)
{
   for (int i = 0; i < height; i++)
   {
      free(p[i]);
   }
   free(p);
}


/*******************************************************************************
 * allocate2Di - allocate 2D array of integers
 *
 * INPUTS: int height - number of rows
 *         int width - number of columns
 *
 * RETURNS: int** - pointer to allocated memory
 *
 * ASSUMES:
 *
 * EFFECTS:  allocates 2D array
 *
 * ERROR CONDITIONS: will return NULL if unable to allocate memory
 *
 * COMMENTS:
 *
 ******************************************************************************/   
int **allocate2Di(int height, int width)
{

   int **p = malloc(sizeof(int*) * height);
   
   if (p == NULL)
      return NULL;
   
   for (int i = 0; i < height; i++)
   {
      p[i] = malloc(sizeof(int) * width);
      if (p[i] == NULL)
      {
         for (int j = 0; j < i; j++)
         {
            free(p[j]);
         }
         return NULL;
      }
   }
   
   return p;

}

/*******************************************************************************
 * free2Di - fee 2D integer array
 *
 * INPUTS: int **iMatrix - 2D memory to free
 *         int height - number of rows
 *
 * RETURNS: void
 *
 * ASSUMES:
 *
 * EFFECTS: frees 2D array
 *
 * ERROR CONDITIONS:
 *
 * COMMENTS:
 *
 ******************************************************************************/
void free2Di(int **p, int height)
{
   for (int i = 0; i < height; i++)
   {
      free(p[i]);
   }
   free(p);
}


/*******************************************************************************
 * fprint_chr - takes chromosome number and prints out printer friendly version
 *
 * INPUTS:   FILE *stream - stream to write to
 *
 * RETURNS:  passes along return value from fprintf
 *
 * ASSUMES:
 *
 * EFFECTS:  characters are written to stream
 *
 * ERROR CONDITIONS: if chr is invalid nothing is written, EHMM is returned
 *
 * COMMENTS:
 *
 ******************************************************************************/
int fprint_chr(FILE *stream, int chr)
{

   if (chr < 1 || chr > CHR_MAX)
   {
      return EHMM;
   }
   else if (chr == CHR_X)
   {
      return fprintf(stream, "X");
   }
   else if (chr == CHR_Y)
   {
      return fprintf(stream, "Y");
   }
   else if (chr == CHR_M)
   {
      return fprintf(stream, "M");
   }
   else
   {
      return fprintf(stream, "%d", chr);
   }


}

/*******************************************************************************
 * expandStrainList - takes comma delimited list that may included ranges
 *    e.g. 1-5) and expands all ranges into comma delimited numbers
 *    1,2,5-10 -> 1,2,5,6,7,8,9,10
 *
 * INPUTS:   char *list - pointer to string containing original list
 *           int num_strains_total - total number of strains
 *           char **expanded_list = pointer to string pointer to hold final 
 *              list
 *           int *num_strains_expanded - pointer to store number of strains in the 
 *               expanded list
 *
 * RETURNS:  0 for sucess, 1 for false
 *
 * ASSUMES:
 *
 * EFFECTS:
 *
 * ERROR CONDITIONS: invalid number in list
 *
 * COMMENTS:
 *
 ******************************************************************************/
int expandStrainList(char *list, int num_strains_total, 
                     char **expanded_list, int *num_strains_expanded)
{
   *expanded_list = malloc(LIST_SIZE);
   char tmp_str[strlen(list) + 1];
   char start_str[5], end_str[5];
   
   *expanded_list[0] = 0;
   start_str[4] = 0;
   end_str[4] = 0;
   
   if (strcmp(list, "all") == 0)
   {
      
      for(int i = 1; i <= num_strains_total; i++)
      {
         sprintf(&((*expanded_list)[strlen(*expanded_list)]), "%d,", i);
         (*num_strains_expanded)++;
      }
         (*expanded_list)[strlen(*expanded_list) - 1] = 0;
         return 0;
   }
   
   /* strtok is destructive so we make a copy */
   strcpy(tmp_str, list);
   
   /*pointers for strtok*/
   char *token, *idx;
   idx = tmp_str;
   
   for (token = strtok(idx, ","); token != NULL;
         token = strtok(NULL, ","))
   {        
      
      
      start_str[0] = (char)0;
      start_str[1] = (char)0;
      end_str[0] = '0';
      end_str[1] = (char)0;

      /* pop characters off until we hit a '-' char or end of string 
         We must take care that we don't over flow our short 4 character
         strings.  (hence the i < 5, i < 4, and j < 4 tests)
         
         valid tokens are one to four digits, followed by an optional 
         dash and one to four more digits
         
         */
         
      int i; /* we need the value of this after the for loop terminates */
 
      /* look at up to the first 5 characters */
      for (i = 0; i < strlen(token) && i < 5; i++)
      {
         /* grab up to four non "-" characters, the 5th spot in the token 
            is reserved for '-' */
         if (token[i] != '-' && i < 4)
         {
            start_str[i] = token[i];
         }
         else if (token[i] == '-')
         {
            i++;
            break;
         }
         else
         {
            /* invalid format - 5th character in token isn't '-' */
            strncpy(*expanded_list, token, LIST_SIZE-1);
            *expanded_list[LIST_SIZE-1] = '\0'; 
            return 1;
         }
      }
      
      for(int j = 0; i < strlen(token) && j < 4; i++,j++)
      {
         end_str[j] = token[i];
      }
     

      int start = (int)strtol(start_str, NULL, 10);
      int end = (int)strtol(end_str, NULL, 10);

      if (end != 0 && start != 0)  /* token was a range of numbers */
      {
         /* validate range */
         if ( !(start < end) || (start < 1) || 
               (start > num_strains_total) || 
               (end > num_strains_total) )
         {
            strncpy(*expanded_list, token, LIST_SIZE-1);
            *expanded_list[LIST_SIZE-1] = '\0'; 
            return 1;
         }
         else
         {
            for(i = start; i<=end; i++)
            {
               snprintf(&((*expanded_list)[strlen(*expanded_list)]), 
                        LIST_SIZE - strlen(*expanded_list), "%d,", i);
               (*num_strains_expanded)++;
            }
         }
      }
      else if ( start != 0 )
      {
         if ((start < 1) || (start > num_strains_total))
         {
            strncpy(*expanded_list, token, LIST_SIZE-1);
            *expanded_list[LIST_SIZE-1] = '\0'; 
            return 1;
         }
         else
         {
            snprintf(&((*expanded_list)[strlen(*expanded_list)]),
                    LIST_SIZE - strlen(*expanded_list), "%d,", start);
            (*num_strains_expanded)++;
         }
      }
      else
      {
         strncpy(*expanded_list, token, LIST_SIZE-1);
         *expanded_list[LIST_SIZE-1] = '\0';
         return 1;
      }
         
   }
   /* trim off the trailing ',' character */
   (*expanded_list)[strlen(*expanded_list)-1] = 0;
   return 0;
}
