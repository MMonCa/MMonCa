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
 *  direct.h
 *  feliks
 *
 *  Created by Ignacio Romero on 10/3/06.
 *
 */

#ifndef _direct_h
#define _direct_h

#include "Math/LinearSolvers/linearsolver.h"
#include "Math/Sparse/compressed.h"

/*
#include "dsp_defs.h"
#include "util.h"
*/

#include <superlu/slu_ddefs.h>
#include <superlu/slu_util.h>
#include <superlu/supermatrix.h>


#include <iostream>

class linearSOE;
class commandLine;
class model;


// superLU solver
class direct : public linearSolver
{
	
public:
	
	direct();
	direct(const commandLine &cl);
	~direct();
	
	void	info(std::ostream &of=std::cout);
	void	initialize(const model& theModel);
	bool	solve(linearSOE &theSOE, double &mflops);


private:
	SuperMatrix *L;
    SuperMatrix *U;
	bool        isPermuted;
	CSCMatrix	A;
	longvector	b;
	
	int    *colPerm;	// an array of ndof ints, the column permutations for SuperLU.
	int    *rowPerm;	// an array of ndof ints, the row permutations for SuperLU.
	int    *elTree;		// a int* of ndof, the elimination tree etree for SuperLU.
	
	
	linearSOE	getTheLinearSOE() const;
	bool		solveReal(linearSOE &theSOE, double &mflops);

};



// superLU multithread solver
class mtdirect : public linearSolver{
	
private:
	SuperMatrix *L;
    SuperMatrix *U;
	bool        isPermuted;
		
	
	linearSOE	getTheLinearSOE() const;
	bool  solveReal(linearSOE &theSOE, double &mflops);

public:
	
	mtdirect();
	mtdirect(const commandLine &cl);
	~mtdirect();
	
	void	info(std::ostream &of=std::cout);
	void	initialize(const model& theModel){}
	bool	solve(linearSOE &theSOE, double &mflops);
};



// sparse cholesky solver for symmetric, positive definite matrices
class cholesky : public linearSolver{
	
private:
	
	double *xp;
	double *bp;
	void   *L;
	int    *perm;		// permutation matrix
	int    *invperm;
	bool    isPermuted;
	SCSCMatrix A;
	longvector b;
	
	linearSOE	getTheLinearSOE() const;
	bool  solveReal(linearSOE &theSOE, double &mflops);
	
public:
	
	cholesky(const commandLine &cl);
	cholesky();
	~cholesky();
	
	void	info(std::ostream &of=std::cout);
	void	initialize(const model& theModel);
	bool	solve(linearSOE &theSOE, double &mflops);
};




#endif
