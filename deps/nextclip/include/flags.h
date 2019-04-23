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
 
#ifndef FLAGS_H_
#define FLAGS_H_
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <global.h>

//Global flags
#define  ALL_OFF  		   	 0

#ifndef SHORT_FLAGS
#define SHORT_FLAGS
#endif

typedef uint16_t Flags;
#define  ASSIGNED			  (1 << 0)  //x0001
#define  VISITED			  (1 << 1)  //x0002 
#define  PRUNED				  (1 << 2)  //x0003 
#define  BRANCH_NODE_FORWARD  (1 << 3)  //x0004 
#define  BRANCH_NODE_REVERSE  (1 << 4)  //x0010
#define  X_NODE			 	  (1 << 5)  //x0020
#define  END_NODE_FORWARD  	  (1 << 6)  //x0040
#define  END_NODE_REVERSE  	  (1 << 7)  //x0080
#define  STARTING_FORWARD	  (1 << 8)  //x0100 Used to mark where a new SOLiD read starts. 
#define  TIP_START            (1 << 9)  //x0200 To mark where the tip clip starts
#define  IGNORE_START_NODE	  (1 << 10) //x0400
#define  CURRENT_PATH_FORWARD (1 << 11) //x0800
#define  CURRENT_PATH_REVERSE (1 << 12) //x1000
#define	 Y_START			  (1 << 13) //x2000 
#define  VISITED_FORWARD  	  (1 << 14) //x4000
#define  VISITED_REVERSE  	  (1 << 15) //x8000

boolean flags_check_for_flag(Flags f, Flags * db);
boolean flags_check_for_any_flag(Flags f, Flags * db);
void flags_action_clear_flags(Flags * db);
void flags_action_set_flag(Flags f, Flags * db);
void flags_action_unset_flag(Flags f, Flags * db);

#endif /* FLAGS_H_ */
