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
 *  stepsolver.h
 *  feliks
 *
 *  Created by Ignacio Romero on Wed Nov 12 2003.
 *  Copyright (c) 2003 Universidad Politecnica de Madrid. All rights reserved.
 *
 *  data definition for nonlinear (and linear) solvers. No functions defined,
 *  all the functions are defined in the corresponding solver (Newton-Raphson,
 *  explicit, modal ...)
 */


#ifndef _stepsolver_h
#define _stepsolver_h
	


#include <iostream>
#include <list>
#include "Math/vector.h"


class assembler;
class commandLine;
class energies;
class integrator;
class linearSOE;
class linearSolver;
class loading;
class model;
class postprocessor;
class tscontrol;



class stepsolver
{
	
public:
	
                        stepsolver();
                        stepsolver(const commandLine& cl);
	virtual             ~stepsolver();
	
	static  void        scan(commandLine &cl, stepsolver **ss);
	virtual void        setLinks(model &m, 
                                 integrator &i, linearSOE &soe, linearSolver &ls, 
                                 tscontrol &tsc, energies &e, postprocessor &pp);
	
	
	
	// Information about the step solver
	assembler&          getAssembler(){return *_theAssembler;}
	const assembler&    getAssembler() const {return *_theAssembler;}
	virtual void        info(std::ostream &of=std::cout);
	void                infoSummary(std::ostream &of=std::cout);
	virtual void        convergenceInfo() const;
	inline int          getNRtrials()				{return NRtrials;}
	inline unsigned     getSteps()					{return nsteps;}
	inline bool         isImplicit()				{return implicit;}
	inline bool         maxNRtrialsExceeded()		{return NRtrials > maxNRtrials;}
	
	// request actions from the step solver
	virtual void        advanceTime(double finaltime, std::ostream &of=std::cout);
	virtual double      adaptiveTimeStep(double finaltime, std::ostream &of=std::cout)=0;
	virtual bool        check(std::ostream &of=std::cout);
	void                imposeConstrainedDofs();
	inline  void        incrementCpuTime(double inc) { cputime += inc;}
	virtual void        incrementStep();
	virtual bool        solveGlobalSystem(std::ostream &of=std::cout);
	virtual bool        solveIteration(std::ostream &of=std::cout);
	virtual bool        solveStep(std::ostream &of=std::cout)=0;
	
	// io functions
	virtual void        messageWhenAdvancingStep(std::ostream &of=std::cout);
    virtual void        messageAfterRemeshing(std::ostream &of=std::cout);
	virtual void        report(std::ostream &of=std::cout);
	virtual void        solutionMessage(std::ostream &of=std::cout, double mflops=0)=0;
		
    
    
protected:
	std::string         name;
	bool                implicit;
	unsigned            nsteps;                     // total number of steps solved
	int                 stepIters;					// total num of iterations in current step
	double              hmindt;						// min value of dt employed in solution
	double              hmaxdt;						// max value of dt employed in solution
    
	double              Eerror;
	double              Rerror;
	double              firstEerror;
	double              firstResid;
	double              cputime;					// time spend in the solution of the analysis
	
	int                 NRtrials;					// NR restarts counter
	int                 maxNRtrials;				// times the NR can be restarted before error
	bool                variableStepSize;
	longvector          tmpVector;					// a temporary vector used in the analysis
	bool                _linesearch;
    
    
	
	model*                  theModel;
	integrator*             theIntegrator;
	linearSOE*              theLinearSOE;
	linearSolver*           theLinearSolver;
	tscontrol*              theTSControl;
	energies*               theEnergy;
	postprocessor*          stPostprocessor;
	assembler*              _theAssembler;    
	
	friend class fractionalStepSolver;
	

    
private:
	stepsolver(const stepsolver& ss);
	
	static void linesearchMessage(const double s, const double Gs, std::ostream& of);
	bool        lineSearch(const longvector& initialResidual, longvector& theSolution);

 };

	
#endif

