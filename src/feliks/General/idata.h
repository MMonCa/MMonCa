/* Finite Element Method Module
 *
 * Author: ignacio.romero@imdea.org
 *
 * Copyright 2014 IMDEA Materials Institute, Getafe, Madrid, Spain,
 *      and       Technical University of Madrid (UPM), Madrid, Spain
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */ 
/*
 * idata.h
 *
 * ignacio romero
 * june 2000
 *
 * idatas are collections of integers that include a counter with the number of integers
 * in the list. When accesssing the individual integers, the list starts from the position 0.
 *
 */

#ifndef _idata_h
#define _idata_h

#ifdef __cplusplus
extern "C" {
#endif
	

#include <stdio.h>
 
/*----------------------------------------------------------------------------------------*/
/*                                   data structure                                       */
/*----------------------------------------------------------------------------------------*/
typedef struct{
  int   maxdata;		/* maximum number of integers that can go in the list */
  int   ndata;			/* number of integers that have been added to the list */
  int  *data;			/* the actual list of integers */
}idataT;

typedef idataT *idata;

/*----------------------------------------------------------------------------------------*/
/*                                     functions                                          */
/*----------------------------------------------------------------------------------------*/

/* create and erase new idata */
idata NewEmptyIdata(void);
idata NewIdata(int length);
void  FillIdata(idata idt, const int length, int* dd);
void  FreeIdata(idata idt);
idata CopyIdata(idata idt);

/* Add data and fill in individual integers in the idata list */
void  AppendToIdata(idata idt, int k);
void  IdataAppend(idata id, int nargs, ...);
void  PutInIdata(idata idt, int theint, int position);

/* Operations on idata */
void IdataEmpty(idata idt);
int  IdataOrder(idata idt);
void IdataAppendIdata(idata tothis, idata appendthis);
int  IdataRemoveRepeated(idata idt);

/* retrive individual integers from the list */
int   GetIntInIdata(int n, idata idt);
int   IdataMaximum(idata idt);
int   IdataMinimum(idata idt);

/* obtain information about the idata */
int   IsIntegerInIdata(idata idt, int i);
int   IdataSize(idata idt);
int   PositionInIdata(idata idt, int i);
int   IdataCompare(idata idt1, idata idt2);

/* Display information */
void PrintIdata(idata idt);
void PrintIdataInFile(FILE *fp, idata idt);
void PrintIdataInBinaryFile(FILE *fp, idata idt);

#ifdef __cplusplus
}
#endif
	

#endif

