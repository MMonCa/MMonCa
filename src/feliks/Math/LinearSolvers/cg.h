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
 *  cg.h
 *  feliks
 *
 *  Created by Ignacio Romero on 10/3/06.
 *  Copyright 2006 Universidad Politécnica de Madrid. All rights reserved.
 *
 */

#ifndef _cg_h
#define _cg_h

#include "Math/LinearSolvers/linearsolver.h"
#include "Math/Sparse/compressed.h"
#include "Math/vector.h"

#include <iostream>
#include <string>

class commandLine;
class linearSOE;
class longvector;	

class cg : public linearSolver
	{
		
	public:
		
		cg(const commandLine &cl);
		cg();
		~cg();
		
		void		info(std::ostream &of=std::cout);
		void		initialize(const model& theModel);
		bool		solve(linearSOE &theSOE, double &mflops);
		
		
	private:
		
		std::string preconditioner;
		double droptol;
		int    *perm;		// permutation matrix
		int    *invperm;
		void   *prec;
		double *bp;
		bool    isPermuted;
		longvector p, r, v, x, z;
		longvector diag;
		
		SCSCMatrix	A;
		longvector	b;
		
		
		bool  solveJacobi(linearSOE &theSOE, double &mflops);
		bool  solveILU(linearSOE &theSOE, double &mflops);
	};



#endif
