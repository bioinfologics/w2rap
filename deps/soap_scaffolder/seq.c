/*
 * seq.c
 *
 * Copyright (c) 2008-2012 BGI-Shenzhen <soap at genomics dot org dot cn>.
 *
 * This file is part of SOAPdenovo.
 *
 * SOAPdenovo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SOAPdenovo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SOAPdenovo.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "stdinc.h"
#include "newhash.h"
#include "kmerhash.h"
#include "extfunc.h"
#include "extvab.h"

/*************************************************
Function:
    writeChar2tightString
Description:
    Writes base in specified position of sequence.
Input:
    1. nt:          the base
    2. tightSeq:        the sequence
    3. pos:         the position
Output:
    None.
Return:
    None.
*************************************************/
void writeChar2tightString ( char nt, char * tightSeq, int pos )
{
	char * byte = tightSeq + pos / 4;

	switch ( pos % 4 )
	{
		case 0:
			*byte &= 63;
			*byte += nt << 6;
			return;
		case 1:
			*byte &= 207;
			*byte += nt << 4;
			return;
		case 2:
			*byte &= 243;
			*byte += nt << 2;
			return;
		case 3:
			*byte &= 252;
			*byte += nt;
			return;
	}
}

/*************************************************
Function:
    getCharInTightString
Description:
    Gets the base in sipcified pos of sequence.
Input:
    1. tightSeq:        the sequence
    2. pos:         the position
Output:
    None.
Return:
    The target base.
*************************************************/
char getCharInTightString ( char * tightSeq, int pos )
{
	char * byte = tightSeq + pos / 4;

	switch ( pos % 4 )
	{
		case 3:
			return ( *byte & 3 );
		case 2:
			return ( *byte & 12 ) >> 2;
		case 1:
			return ( *byte & 48 ) >> 4;
		case 0:
			return ( *byte & 192 ) >> 6;
	}

	return 0;
}

/*************************************************
Function:
    reverseComplementSeq
Description:
    Gets the reverse complement of  a sequence.
Input:
    1. seq:         the sequence
    2. len:         the length of the sequence
Output:
    1. bal_seq:     the reversed complement of the sequence
Return:
    None.
*************************************************/
void reverseComplementSeq ( char * seq, int len, char * bal_seq )
{
	int i, index = 0;

	if ( len < 1 )
	{
		return;
	}

	for ( i = len - 1; i >= 0; i-- )
	{
		bal_seq[index++] = int_comp ( seq[i] );
	}

	return;
}
