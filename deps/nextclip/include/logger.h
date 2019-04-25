/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo 
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
 * Author:  Richard Leggett                                             *
 *          The Sainsbury Laboratory, Norwich, UK                       *
 *          richard.leggett@tsl.ac.uk                                   *
 * History: 19-Oct-10: RML: Pulled together from ad hoc code!           *
 *           6-Feb-12: RHRG:Added a function to print to any FILE       *
 *----------------------------------------------------------------------*/

int log_start(char* filename);
void log_printf(char* fmt, ...);
void log_and_screen_printf(char* fmt, ...);
void log_newline(void);
void log_get_timestamp(char* buffer);
void log_write_timestamp(int newline);
void log_progress_bar(int percentage);
void log_set_screen(FILE * f);
