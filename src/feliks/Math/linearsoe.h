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
/* linearsoe.h
 *
 * linear system of equations of the form Ax=b.
 * This objects wrap the matrix, the rhs, the solution, and some extra data, but does not own them
 *
 * ignacio romero
 * june 2000
 */


#ifndef _linearsoe_h
#define _linearsoe_h
		

#include <string>
#include <cstdio>
#include <iostream>

#include "Math/Sparse/sparse.h"
#include "Math/vector.h"

class model;
	


class linearSOE{
	
	public:
	
	// functions to create new system of equations, including an empty matrix A and rhs
	linearSOE();
	linearSOE(model &m, const storageT st, const std::string& dumpfmt="");
	linearSOE(sparseMatrix &A, longvector &b);
	~linearSOE();
	void setLinks(sparseMatrix &Ac, longvector& bc);
	
	// functions to return either the matrix or the rhs
	inline sparseMatrix& getA() {return *A;}
	inline longvector&	 getB() {return *b;}


	// functions to empty the matrix and the rhs
	void  zeroA();
	void  zeroB();


	// functions to get information about the SOE
	inline int    getDimension() {return dimension;}
	void          info(std::ostream &of=cout);


private:
	int           dimension;
	sparseMatrix  *A;			// a pointer to a matrix elsewhere, never allocated
	longvector	  *b;			// a pointer to a vector elsewhere, never allocated
	int           *Rprec;		// row permutation matrix
	int           *Cprec;		// column permutation matrix
	bool			deleteMembersWhenFinished;
	
	friend class petscsolver;
};
	

#endif
