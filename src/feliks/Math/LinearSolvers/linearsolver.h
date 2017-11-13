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
 *  linearSolver.h
 *  ris
 *
 *  Created by Ignacio Romero on Wed May 05 2004.
 *  Copyright (c) 2004 Universidad Politecnica de Madrid. All rights reserved.
 *
 */


#ifndef _linearSolver_h
#define _linearSolver_h


#include <string>
#include "Math/linearsoe.h"


class commandLine;
class model;


class linearSolver
{
	
public:
	
	linearSolver();
	linearSolver(const commandLine& cl);
	linearSolver(const linearSolver& ls);
	virtual ~linearSolver();
	
	
	static void			scan(const commandLine &cl, linearSolver **ls);
	inline std::string	getName() const;
	inline  storageT	getStorage() const;
	linearSOE&			getTheLinearSOE() {return theLinearSOE;};
	virtual void		initialize(const model& theModel)=0;
	virtual void		info(std::ostream &of=cout);
	
	bool				solve(double &cputime, double &mflops);
	bool				solve(double &mflops);
	bool				solve(linearSOE &theSOE, double &cputime, double &mflops);
	virtual bool		solve(linearSOE &theSOE, double &mflops)=0;
	inline  void		setRelTolerance(const double rt);
	
	
protected:
	std::string	name;
	double		rtol;			// relative for convergence
	int         nprocessors;
	storageT    storage;		// this is ugly, but for the moment, solvers need a specific type of sparse matrix
	linearSOE	theLinearSOE;
};



// inlines
inline std::string	linearSolver :: getName() const						{return name;}
inline storageT		linearSolver :: getStorage() const					{return storage;}
inline void			linearSolver :: setRelTolerance(const double rt)	{rtol=rt;}


#endif
