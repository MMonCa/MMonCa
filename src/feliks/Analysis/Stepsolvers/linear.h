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
 *  linearss.h
 *  feliks
 *
 *  Created by Ignacio Romero on dec 22 2005
 *  Copyright (c) 2003 Universidad Politecnica de Madrid. All rights reserved.
 *
 */


#ifndef _linearss_h
#define _linearss_h
	
#include "Analysis/Stepsolvers/stepsolver.h"
#include <iostream>

class commandLine;

class linearss : public stepsolver {

private:
	
	double   Etolerance;					// convergence tolerance in energy norm
	double   Rtolerance;					// convergence tolerance in residual norm	
	

	
public :
	linearss(commandLine &cl);
	linearss();
	~linearss();
	
	bool          solveStep(std::ostream &of=std::cout);
	inline double getFirstResidual()		{return firstResid;}
	inline double getResidualError()		{return Rerror;}
	virtual void	info(std::ostream &of=std::cout);
	
	// request actions from the step solver
	double adaptiveTimeStep(double finalTime, std::ostream &of=std::cout);
	void   solutionMessage(std::ostream &of=std::cout, double mflops=0.0);
 };


#endif
