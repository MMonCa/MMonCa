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
 *  newtonraphson.h
 *  feliks
 *
 *  Created by Ignacio Romero on dec 22 2005
 *  Copyright (c) 2003 Universidad Politecnica de Madrid. All rights reserved.
 *
 */


#ifndef _newtonraphson_h
#define _newtonraphson_h
	
#include "Analysis/Stepsolvers/stepsolver.h"
#include "Io/usercommand.h"
#include <iostream>


class newtonraphson : public stepsolver {

private:
	
	int     iteration;					// num of nr iterations in current nr step
	int     totalIterations;			// total number of iterations in solution history
	int     maxIterations;				// max number of iterations allowed in 1 t-step
	double  ErelTolerance;				// convergence tolerance in energy norm
	double  EabsTolerance;
	double  Rtolerance;					// convergence tolerance in residual norm
	int		recomputeTangent;			// every certain number of iterations
	double  divergenceRerror;
	int		preiterations;
	
	friend class fractionalStepSolver;
	
public :
	newtonraphson(const commandLine &cl);
	newtonraphson();
	~newtonraphson();
	
	virtual void	info(std::ostream &of=std::cout);
	bool			solveStep(std::ostream &of=std::cout);
	inline double	getFirstResidual() const	{return firstResid;}
	inline int		getNRtrials() const		{return NRtrials;}
	inline double	getResidualError() const	{return Rerror;}
	inline int		getStepIterations() const	{return stepIters;}
	virtual void	report(std::ostream &of=std::cout);

	inline double& energyErrorTolerance()   {return ErelTolerance;}
	
	// request actions from the step solver
	double adaptiveTimeStep(double finalTime, std::ostream &of=std::cout);
	void   solutionMessage(std::ostream &of=std::cout, double mflops=0.0);
	bool   solveIteration(std::ostream &of=std::cout);
	
 };


#endif
