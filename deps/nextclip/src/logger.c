/*
 *Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 *  
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * Development team: 
 *       R. Ramirez-Gonzalez (Ricardo.Ramirez-Gonzalez@bbsrc.ac.uk)
 *       R. Leggett (richard@leggettnet.org.uk)
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */

/*----------------------------------------------------------------------*
 * File:    logger.c                                                    *
 * Purpose: Handle debugging/log files.                                 *
 * Authors: Richard Leggett                                             *
 *          The Sainsbury Laboratory, Norwich, UK                       *
 *          richard.leggett@tsl.ac.uk    								*
 * 	        Ricardo H. Ramirez G.										*
 * 			The Genome Analysis Centre                	       			*
 * 			ricardo.ramirez-gonzalez@bbsrc.ac.uk						*		
 * History: 19-Oct-10: RML: Pulled together from ad hoc code!  			*
 *          25-Oct-10: Added progress bar								*
 *----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <stdarg.h>
#include <logger.h>

char* log_filename = NULL;
static FILE * screen = NULL;

/*----------------------------------------------------------------------*
 * Function: log_start                                                  *
 * Purpose:  Create a new log file, ready for writing.                  *
 * Params:   filename -> name of log file.                              *
 * Returns:  1 for success, 0 if file couldn't be opened.               *
 *----------------------------------------------------------------------*/
int log_start(char* filename)
{
	int rc = 1;
	FILE *fp;

	if (!filename) {
		rc = 0;
	} else {	
		log_filename = filename;
		fp = fopen(log_filename, "w");

		if (fp) {
			char timestring[64];
			log_get_timestamp(timestring);
			fprintf(fp, "%s Log started\n", timestring);
			fclose(fp);
		} else {
			log_filename = 0;
			rc = 0;
		}
	}
	
	return rc;
}

/*----------------------------------------------------------------------*
 * Function: log_printf                                                 *
 * Purpose:  Like printf, but outputs sent to log file.                 *
 * Params:   fmt -> printf format string.                               *
 *           + optional parameters.                                     *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void log_printf(char* fmt, ...)
{
    va_list args;

	if (log_filename) {
		FILE *fp = fopen(log_filename, "a");
		if (fp) {
			va_start(args, fmt);
			vfprintf(fp, fmt, args);
            fclose(fp);
		}
	}
}

/*----------------------------------------------------------------------*
 * Function: log_and_screen_printf                                      *
 * Purpose:  Like printf, but outputs sent to log file AND screen.      *
 * Params:   fmt -> printf format string.                               *
 *           + optional parameters.                                     *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void log_and_screen_printf(char* fmt, ...)
{
	va_list args;

	// Output to log
	if (log_filename) {
		FILE *fp = fopen(log_filename, "a");
		if (fp) {
			va_start(args, fmt);
			vfprintf(fp, fmt, args);
            fclose(fp);
		}
	}
	
	// Output to screen
    if(screen == NULL){
        screen = stdout;
    }
	va_start(args, fmt);
	vfprintf(screen, fmt, args);
}

/*----------------------------------------------------------------------*
 * Function: log_newline                                                *
 * Purpose:  Print a newline to the log file.                           *
 * Params:   None.                                                      *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void log_newline(void)
{
	log_printf("\n");
}

/*----------------------------------------------------------------------*
 * Function: log_get_timestamp                                          *
 * Purpose:  Make a timestamp string.                                   *
 * Params:   buffer -> char array to store string.                      *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void log_get_timestamp(char* buffer)
{
	time_t ltime = time(NULL);
	char *timestring = asctime(localtime(&ltime));
	strcpy(buffer, timestring);
	buffer[strlen(timestring)-1] = 0;
}

/*----------------------------------------------------------------------*
 * Function: log_write_timestamp                                        *
 * Purpose:  Write a timestamp to the logfile.                          *
 * Params:   newline = 1 to follow timestamp with a newline character.  *
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/
void log_write_timestamp(int newline)
{
    
	char timestring[64];
	log_get_timestamp(timestring);
    
	if (newline) {
	    log_printf("----- %s -----\n", timestring);
            //log_and_screen_printf("----- %s -----\n", timestring);
	} else {
	    log_printf("%s", timestring);
            //log_and_screen_printf("%s", timestring);
	}
}

/*----------------------------------------------------------------------*
 * Function: log_progress_bar                                			*
 * Purpose:  Print a progress var at the default output.                          	*
 * Params:   percentage.Integer represnting the progress of the process	*
 * Returns:  None.                                                      *
 *----------------------------------------------------------------------*/

void log_progress_bar(int percentage)
{
    if(screen == NULL){
        screen = stdout;
    }
#ifndef SUPRRESS_PROGRESS_BAR
	int i=0;
	fprintf(screen, "\r%i%%\t[",percentage);
	for(i=0;i<100;i++)
		if(i<percentage)
			fprintf(screen, "=");
		else
			fprintf(screen, " ");
	fprintf(screen, "]");
	fflush(screen);
#endif
}

void log_set_screen(FILE * f){
    screen = f;
}
