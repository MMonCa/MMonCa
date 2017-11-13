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
/* longvector.h
 * iro
 * october 1999 
 
 *
 */


#ifndef _longvector_h
#define _longvector_h


#include "General/idata.h"
#include <iostream>
#include <vector>



class longvector
{

    
public:

	int     length;
    double *data;

	longvector(const size_t l);
	longvector();
	longvector(const longvector& v);
	~longvector();

	longvector&         operator=(const longvector& v);
	longvector&         operator+=(const longvector& v);
	longvector&         operator-=(const longvector& v);
	longvector&         operator*=(const double a);
    longvector          operator+(const longvector& v);
    longvector          operator-(const longvector& v);

	inline			double&		operator[](const int i)			{return data[i];}
	inline const	double&		operator[](const int i) const	{return data[i];}
	inline			double&		operator()(const int i)			{return data[i];}
	inline const	double&		operator()(const int i) const	{return data[i];}

	
 	longvector  subVector(const int start, const int end) const;
	void		changeSign();
	void		resize(const int length);
	void		scale(const double s);
	void        wrap(double *data, int length);
	void		invertComponents();
	static void	multiplyComponentwise(const longvector& a, const longvector& b, longvector& ab);
	void		unwrap();
 

	// math operations on longvectors
	//void    CopySubVector(longvector v1, longvector v2, int start, int end); /* put v1 into v2 */
	double		dot(const longvector& v2) const;
	double		norm() const;
	double		squaredNorm() const;


	// functions to get information from a longvector
	void		print(std::ostream& of=std::cout);
	int			getLength() const;


	/* functions to modify longvectors */
	void    assembleInVector(const longvector& block, idata list);
	void    assembleInVector(const longvector& block, const std::vector<int>& list);

	void    fill(double *data , int ndata);
	void    setZero();
};	




#ifdef __cplusplus
extern "C" {
#endif
	

/* operations on longvectors, as arrays of doubles */
void    Dcopy(const double *from, double *to, const int l);
double  Ddot(double *v1, double *v2, const int l);
double  Dnorm(double *v1, const int l);
double  Dnormalize(double *v1, const int l);
void    Dprint(double *v, const int l);
void    Dzero(double *v, const int l);


/* operations on arrays of integers */
void Icopy(const int *from, int *to, const int l);
void Izero(int *v, const int l);
void Iprint(int *v, const int l);


#ifdef __cplusplus
}
#endif


#endif
