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
/* longvector.cpp
 * longvector manipulation functions
 *
 * ignacio romero, oct 99
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <cassert>

#include <math.h>
#include "Math/vector.h"
#include "General/idata.h"
#include "Main/feliks.h"

#define NUM_PER_LINE  5          /* numbers per line, for screen output */


longvector :: longvector() :
	length(0),
	data(0)
{
}




longvector :: longvector(const size_t l) :
	length((int)l),
	data(0)
{
	data = new double[l];
	Dzero(data, (int)l);
}


longvector :: longvector(const longvector& v) :
	length( v.length ),
	data(0)
{
	if (length > 0)
	{
		data = new double[length];
        memcpy(data, v.data, sizeof(double)*length);
	}
}



longvector :: ~longvector()
{
	if (data != 0) delete [] data;
}



/* assembles the data in the longvector 'block' inside the longvector 'big' in
the positions indicated by 'list'. An entry '-1' in list indicates
that the corresponding data need not be assembled.
*/
void  longvector :: assembleInVector(const longvector& block, idata list)
{
	int m, globalSlot, len;
	
	double *Big   = data;
	double *Block = block.data;
	
	len = list->ndata;
	for (m=0; m<len; m++)
    {
		globalSlot = list->data[m];
		if(globalSlot >=0) Big[globalSlot] += Block[m];
    }
}


void  longvector :: assembleInVector(const longvector& block, const std::vector<int>& list)
{
	int m, globalSlot, len;
	
	double *Big   = data;
	double *Block = block.data;
	
	len = list.size();
	for (m=0; m<len; m++)
    {
		globalSlot = list[m];
		if(globalSlot >=0) Big[globalSlot] += Block[m];
    }
}




void longvector ::	changeSign()
{
	for (int i=0; i<length; i++)
		data[i] = -data[i];
}



double longvector :: dot(const longvector& v2) const
{
    if ( length != v2.length )
	{
        printf("\n Error in Dot().");
        printf("\n Different size longvectors");
	}
	
    double dot=0.0;
	for (int k=0; k<length; k++) 
        dot += data[k]*v2.data[k];

	return dot;
}




void  longvector :: fill(double *cdata , int ndata)
{
	this->resize(ndata);
    memcpy(this->data, cdata, sizeof(double)*ndata);
}



void longvector :: invertComponents()
{
	for (int k=0; k<length; k++) data[k] = 1.0/data[k];
}



longvector& longvector ::  operator=(const longvector &v)
{
	this->resize(v.length);
    memcpy(this->data, v.data, sizeof(double)*v.length);
	return *this;
}


longvector&	longvector :: operator+=(const longvector &v)
{
	for (int i=0; i< length; i++) 	data[i] += v.data[i];
	return *this;
}


longvector longvector :: operator+(const longvector& v)
{
    assert( v.length == this->length);
    longvector sum(length);
	for (int i=0; i< length; i++)
        sum[i] = (*this)[i] + v[i];
    
    return sum;
}



longvector longvector :: operator-(const longvector& v)
{
    assert( v.length == this->length);
    longvector dif(length);
	for (int i=0; i< length; i++)
        dif[i] = (*this)[i] - v[i];
    
    return dif;
}



longvector&	longvector :: operator-=(const longvector &v)
{
	for (int i=0; i< length; i++) 	data[i] -= v.data[i];
	return *this;
}



longvector&	longvector :: operator*=(const double a)
{
	for (int i=0; i< length; i++) 	data[i] *= a;
	return *this;
}


void longvector :: print(std::ostream& of)
{
	int i, pstart, pend;
	
	if (length > 0 )
    {
		of << "\n" << std::setprecision(6);
		pstart = 0;
		pend   = std::min<int>(length, NUM_PER_LINE);
		while (pstart < length)
		{
			of << "\n" <<  std::setw(3) << pstart+1 << "-" << std::setw(3) << pend;
			
			// column loop
			for (i=pstart; i<pend; i++) of << std::right << std::scientific << std::setw(of.precision()+8) << data[i];			
			
			pstart += NUM_PER_LINE;
			pend    = std::min<int>(length , pend + NUM_PER_LINE);
		}
		of << std::flush;
    }
}



void longvector :: resize(const int newlength)
{
	if ( this->length != newlength)
    {
		delete [] data;
		this->data = new double[newlength];
    	this->length = newlength;
	}
}




void longvector :: scale(const double s)
{
	for (int i=0; i<length; i++)	data[i] *= s;
}



longvector longvector :: subVector(const int start, const int end) const
{
  longvector subvector;
    
  if ( start<0 || end>length)
      printf("\n ERROR in Sublongvector. Can not extract such a longvector");

  else
    {
      subvector.resize(end-start+1);
      for (int i=start; i<=end; i++)
		  subvector.data[i-start] = data[i];
    }

  return subvector;
}




int longvector :: getLength() const
{
  return length;
}


void longvector :: multiplyComponentwise(const longvector& a, const longvector& b, longvector& ab)
{
	for (int i=0; i<a.length; i++)	
		ab.data[i] = a.data[i] * b.data[i];
}




double longvector :: norm() const
{
	return sqrt( squaredNorm() );
}



double longvector :: squaredNorm() const
{
	return this->dot(*this);
}


void longvector :: unwrap()
{
	length = 0;
	data   = 0;
}


/* this function takes a double *, which should not be null and wraps it
   inside a longvector structure, indicating its length. This allows to 
   change the length of a longvector without changing its allocated data and
   to have longvectors which are statically allocated, avoiding having to
   allocate large chunks of memory every time one of them is used 
*/
void   longvector :: wrap(double *data_, int length_)
{
  data     = data_;
  length   = length_;
}




/* initialized to zero a longvector */
void longvector :: setZero()
{
	Dzero(data, length);
}


/*
void Zcopy(complex double *from, complex double *to, const int l)
{
	Dcopy((double *) from, (double *) to, 2*l);
}
*/
/*
double Zdot(complex *v1, complex *v2, const int l)
{
	complex  d=0+0*I;
	int    i;
	
	for (i=0; i<l; i++) d += v1[i]*conj(v2[i]);
	
	return norm(d);
}


void Zzero(complex double *v, const int l)
{
	Dzero((double *)v, 2*l);
}
*/





void Dcopy(const double *from, double *to, const int l)
{
	for (int i=0; i<l; i++) to[i] = from[i];
}





double Ddot(double *v1, double *v2, const int l)
{
	double d=0.0;
	for (int k=0; k<l; k++) d += v1[k]*v2[k];
	return d;
}





/* compute the 2-norm of a longvector */
double  Dnorm(double *v, const int l)
{
	return sqrt(Ddot(v, v, l));
}


/* divide a longvector by its 2-norm */
double  Dnormalize(double *v, const int l)
{
    double n, in;
    int    i;
	extern double global_macheps;
    
    n  = Dnorm(v, l);
	if ( n > 2.0*l*global_macheps)
	{
		in = 1.0/n;
		for (i=0; i<l; i++) v[i] *= in;
	}
	else
		Dzero(v, l);
	
    return n;
}


void    Dprint(double *v, const int l)
{
	longvector tmp(l);
	tmp.fill(v, l);
	tmp.print();
}


void Dzero(double *v, const int l)
{
    for (int k=0; k<l; k++) v[k] = 0.0;
}




void Icopy(const int *from, int *to, const int l)
{
	int i;
	for (i=0; i<l; i++) to[i] = from[i];
}


void Iprint(int *v, const int l)
{
	idataT tmp;
	tmp.ndata = tmp.maxdata = l;
	tmp.data  = v;
	PrintIdata(&tmp);
}


void Izero(int *v, const int l)
{
	int i;
	for (i=0; i<l; i++) v[i] = 0;
}

/*
int InsertVectorInVector(longvector& v1,longvector& v2,int position)
{
	int  len1, len2;
	int  ret, i;
	
	
	len1 = v1.length;
	len2 = v2.length;
	
	if ( position+len1 > len2)
    {
		printf("\n Error in InsertVectorInVector()");
		printf("\n You can not insert a longvector so large.");
		ret  = -1;
    }
	else
    {
		for (i=position; i<position+len1; i++)
			*(v2.data+i) = *(v1.data+(i-position));
		ret = 1;
    }
	return(ret);
}

*/

