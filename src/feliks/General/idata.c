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
 * idata.c
 *
 * ignacio romero
 * june 2000
 *
 * idatas are collections of integers
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "General/idata.h"
#include "Io/message.h"


#define IDATA_SIZE 10		/* default idata length */
#define IRISMAX(a,b)   ( (a>b) ? a : b )
#define IRISMIN(a,b)   ( (a<b) ? a : b )


/* private functions */
static int IdataOrderAux(const void *int1, const void *int2);



void IdataAppendIdata(idata tothis, idata appendthis)
{
    int k;
    
    for (k=0; k< appendthis->ndata; k++) AppendToIdata(tothis, appendthis->data[k]);
}



void AppendToIdata(idata id, int i)
{
    if (id == NULL) ErrorMessage("in AppendIntToIdata. NULL idata");
    
    // stretch idata if necessary
    if (id->ndata == id->maxdata)
    {
		int *newspace=NULL, k;
		int newsize=0;
		
        /* create more space and copy old stuff in it */
        newsize   = id->ndata + (id->ndata/2) + 1;
        newspace  = (int*) malloc(sizeof(int)*newsize);
        for (k=0; k<id->ndata; k++) newspace[k] = id->data[k];
        
        /* free the old space and link to larger idata*/
        free(id->data);
        id->data    = newspace;
        id->maxdata = newsize;
    }
    
    id->data[id->ndata] = i;
    (id->ndata)++;
    
}


void IdataAppend(idata id, int nargs, ...)
{
    va_list args;
    int k, v;
    
    if (id == NULL) ErrorMessage("in AppendIntToIdata. NULL idata");
    
    //call with the last named argument
    va_start(args, nargs);
    
    for (k=0; k<nargs; k++)
    {
        v = va_arg(args, int);
        AppendToIdata(id, v);
    }
    va_end(args);
}




idata CopyIdata(idata id)
{
  idata id2=NULL;
  int   a;

  if (id == NULL) 
    return NULL;
  else
    {
      id2 = NewIdata(id->maxdata);
      for (a=0; a<id->ndata; a++) id2->data[a] = id->data[a];
      id2->ndata = id->ndata;
   }
  return id2;
}


void  FillIdata(idata idt, const int length, int* dd)
{
	IdataEmpty(idt);
	if ( idt->maxdata < length)
	{
		free(idt->data);
		idt->data   = (int*) malloc( length * sizeof(int) );
		idt->maxdata= length;
	}

	idt->ndata = length;
	
	int i;
	for (i=0; i<idt->ndata; i++) idt->data[i] = dd[i];
}


void FreeIdata(idata id)
{
  free(id->data);
  free(id);
}



// compares two idata, as a dictionary. This means that if two idatas
// have different lengths, but have different sizes, the shorter one
// is "smaller"
// function returns -1, 0, +1
// +1: the first one is larger
//  0: they are identical
// -1: the first one is smaller,
int IdataCompare(idata idt1, idata idt2)
{
   int k, l1, l2, l, ret;
    
    l1 = idt1->ndata;
    l2 = idt2->ndata;
    l  = IRISMIN(l1, l2);
    
    ret = 0;
    for (k=0; k<l; k++)
    {
        if (idt1->data[k] < idt2->data[k])
        {
            ret = -1;
            break;
        }
        else if (idt1->data[k] > idt2->data[k])
        {
            ret = 1;
            break;
        }
    }
    
    if      (ret == 0 && l1<l2)  ret = -1;
    else if (ret == 0 && l1>l2)  ret =  1;
    
    return ret;
}


void IdataEmpty(idata idt)
{
	idt->ndata = 0;
}



int   IdataMaximum(idata idt)
{
	int i, imax;

	imax = idt->data[0];
	for (i=1; i<idt->ndata; i++) imax = IRISMAX(imax, idt->data[i]);
	
	return imax;
}


int  IdataMinimum(idata idt)
{
	int i, imin;
	
	imin = idt->data[0];
	
	for (i=1; i<idt->ndata; i++) imin = IRISMIN(imin, idt->data[i]);
	
	return imin;
}	



/*----------------------------------------------------------------------------------------*/
/*                                      IdataSize                                         */
/*                                                                                        */
/* Purpose:   retrieve the number of integers in the idata list                           */
/*                                                                                        */
/* Input:     the idata structure                                                         */
/* Output:    the number of integers in the idata list                                    */
/*                                                                                        */
/*----------------------------------------------------------------------------------------*/
int IdataSize(idata id)
{
    if (id==NULL) ErrorMessage("in IdataSize. NULL idata");
    return(id->ndata);
}





int  IdataRemoveRepeated(idata idt)
{
    int different, repeated;
    int last, k, *new=NULL;
    
    // quick return
    if (idt->ndata == 0) return 0;
    
    // first we need to order the contents of the list
    IdataOrder(idt);
    
    // now we traverse it checking for repeated values, counting the unique ones
    // count the number of different faces  
    different = 1;
    repeated  = 0;
    last      = idt->data[0];
    for (k=1; k<idt->ndata; k++)
    {
        if (idt->data[k] != last)
        {
            different++;
            last = idt->data[k];
        }
        else repeated++;
    }
    
    // allocate space for non repeated values and fill up the array
    new = (int*) malloc(sizeof(int)*different);
    
    different=0;
    new[different++] = idt->data[0];
    last   = new[0];
    for (k=1; k<idt->ndata; k++)
    {
        if (idt->data[k] != last)
        {
            new[different++] = idt->data[k];
            last = idt->data[k];
        }
    }
    
    // we clean old Idata and replace it with the unique data
    free(idt->data);
    idt->data  = new;
    idt->ndata = different;
    idt->maxdata = different;
    
    return repeated;
}




/*----------------------------------------------------------------------------------------*/
/*                                      GetIntInIdata                                     */
/*                                                                                        */
/* Purpose:   returns the integer that accupies a certain slot in the idata list          */
/*                                                                                        */
/* Input:     n : the position. Starting from 0.                                          */
/*            id: the idata                                                               */
/*                                                                                        */
/* Output:    the integer in idata                                                        */
/*                                                                                        */
/*----------------------------------------------------------------------------------------*/
int GetIntInIdata(int n, idata id)
{
    if (id == NULL) ErrorMessage("in GetIntInIdata. NULL idata");
    if (n  >= id->ndata) ErrorMessage("in GetIntInIdata. Idata is not so large.");

    return( id->data[n] );
}




/* return boolean to indicate if an integer is part of an idata. It can be used
   as in
      if (IsIntegerInIdata(id, i)) ...
 */
int   IsIntegerInIdata(idata id, int i)
{
    int n, a;

    n = id->ndata;

    for (a=0; a<n; a++)
    {
        if (id->data[a] == i) return 1;
    }

    return 0;
}


/*----------------------------------------------------------------------------------------*/
/*                                      NewEmptyIdata                                     */
/*                                                                                        */
/* Purpose:   creates new idata of standard length with no data in it                     */
/*                                                                                        */
/* Output:    the new idata structure                                                     */
/*                                                                                        */
/*----------------------------------------------------------------------------------------*/
idata NewEmptyIdata(void)
{
    idata  id;

    id          = (idata) malloc(sizeof *((idata) NULL));
    id->ndata   = 0;
    id->maxdata = IDATA_SIZE;
    id->data    = (int *) malloc(IDATA_SIZE * sizeof(int));

    return id;
}



idata NewIdata(int n)
{
    idata id;

    id          = (idata) malloc( sizeof *((idata) NULL));
    id->data    = (int *) malloc( n*sizeof (int) );
    id->ndata   = 0;
    id->maxdata = n;

    return id;
}



/* orders the integers in an idata in increasing order, using the qsort
 * implementation from the standard library. For this we need an auxiliary
 * function IdataOrderAux. See below
 */
int  IdataOrder(idata id)
{
    qsort(id->data, id->ndata, sizeof(int), IdataOrderAux);
    return 1;
}


static int IdataOrderAux(const void *int1, const void *int2)
{
    int *i1, *i2, ret;

    i1 = (int *) int1;
    i2 = (int *) int2;

    if (*i1 < *i2)
        ret = -1;
    else if (*i1 == *i2)
        ret = 0;
    else
        ret = 1;

    return ret;
}
    



/* 
 * returns the position of the first appearence of 'i' inside 'id'. If it doesn't appear, 
 * it returns -1
 */
int PositionInIdata(idata id, int i)
{
  int pos, k;

  pos = -1;

  if (id != NULL) 
    {
    for (k=0; k<id->ndata; k++)
       if (id->data[k] == i) 
	  {
	    pos = k;
	    break;
	  }
    }

  return pos;
}



void PrintIdata(idata id)
{
    int k;

    if ( id == NULL) ErrorMessage("in PrintIdata. Null idata.");

    for (k=0; k< id->ndata; k++) printf("    %5d", id->data[k]);
}


void PrintIdataInFile(FILE *fp, idata id)
{
    int k;

    if (id == NULL) ErrorMessage ("in PrintIdataInFile. Null idata.");
    if (fp == NULL) ErrorMessage ("in PrintIdataInFile. Closed file.");

    for (k=0; k< id->ndata; k++) fprintf(fp,"  %5d", id->data[k]);
}


void PrintIdataInBinaryFile(FILE *fp, idata id)
{

    if (id == NULL) ErrorMessage ("in PrintIdataInFile. Null idata.");
    if (fp == NULL) ErrorMessage ("in PrintIdataInFile. Closed file.");

	fwrite(id->data, sizeof(int), id->ndata, fp);
}



void PutInIdata(idata id, int theint, int position)
{
    if (id == NULL) ErrorMessage("in PutInIdata. NULL idata");
    if (position >= id->maxdata) ErrorMessage("in PutInIdata. Idata is not large enough.");

    id->data[position] = theint;
    id->ndata = position+1 > id->ndata ? position+1 : id->ndata;
}



/*        E x a m p l e
 *

       idata id;
           .
           .
           .
       id = NewIdata(3);
       PutInIdata(id, 11, 0);
       PutInIdata(id, 22, 1);
       PutInIdata(id, 33, 2);
       AppendToIdata(id, 99);
       AppendToIdata(id, 99);
       AppendToIdata(id, 99);
       PrintIdata(id);
           .
           .
           .
 */
