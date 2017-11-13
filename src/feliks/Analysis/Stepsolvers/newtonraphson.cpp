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
 *  newtonraphson.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on Thu Sep 11 2003.
 *  Copyright (c) 2003 Universidad Politecnica de Madrid. All rights reserved.
 *
 */

#include "Analysis/Stepsolvers/newtonraphson.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include "boost/foreach.hpp"

#include "Analysis/analysis.h"
#include "Analysis/assembler.h"
#include "Analysis/energies.h"
#include "Main/feliks.h"
#include "Math/linearsoe.h"
#include "Io/logger.h"

#include "Model/model.h"
#include "Model/Parts/modelpart.h"

#include "Io/message.h"
#include "Math/Sparse/sparse.h"


newtonraphson :: newtonraphson() :
	iteration(0),
	maxIterations(MAX_ITERATIONS),
	totalIterations(0),
	recomputeTangent(1),
	ErelTolerance(ENERGY_TOLERANCE),
	EabsTolerance(ENERGY_ABS_TOLERANCE),
	Rtolerance(RESIDUAL_TOLERANCE),
	divergenceRerror(1e10),
	preiterations(0)
{
	implicit    = true;
	name        = "Newton-Raphson";
}


newtonraphson :: newtonraphson(const commandLine &cl) :
	stepsolver(cl),
	iteration(0),
	maxIterations(MAX_ITERATIONS),
	recomputeTangent(1),
	ErelTolerance(ENERGY_TOLERANCE),
	EabsTolerance(ENERGY_ABS_TOLERANCE),
	Rtolerance(RESIDUAL_TOLERANCE),
	divergenceRerror(1e10),
	preiterations(0)
{
	implicit    = true;
	name        = "Newton-Raphson";
	
	for (int k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
		if       ( uc.keyword() == "maxiterations")     maxIterations   = static_cast<int>(uc.value());
		else if  ( uc.keyword() == "maxtrials" )		maxNRtrials     = static_cast<int>(uc.value());
		else if  ( uc.keyword() == "update" )           recomputeTangent= static_cast<int>(uc.value());
		else if  ( uc.keyword() == "etolerance")		ErelTolerance	= uc.value();
		else if  ( uc.keyword() == "abstolerance")      EabsTolerance	= uc.value();
		else if  ( uc.keyword() == "preiterations")     preiterations	= static_cast<int>(uc.value());
	}		
}



newtonraphson :: ~newtonraphson()
{
}




double newtonraphson :: adaptiveTimeStep(double finalTime, std::ostream& of)
{
	double newdt, t = theIntegrator->theTime, dt = theTSControl->dt;
	extern double global_macheps;
	
	theTSControl->advanceDtVariables();
	
	if (t == 0.0 && dt == 0.0)
	{
		newdt = 1000.0/theModel->maxEigenvalueEstimate(); // the 1000 factor is a guess
		Message("\n  No starting time step size. Guessing with %e", newdt);
	}
    else if (t == 0.0 && dt > 0.0)
		newdt = dt;
	else
		newdt = theTSControl->adaptiveForwardDt();
    
	
	// final check. Do not go over final time
	// some extra space to allow for fp errors
	double gap = (finalTime-t) + finalTime*0.5*(nsteps+1)*global_macheps;
	if (finalTime != -1 && newdt-gap > 0.0)
    {
		of << "\n  Correcting NR time step not to exceed final time, dt = " << scientific << newdt << ", gap = " << gap;
        newdt = finalTime - theIntegrator->theTime;
		of << "\n  Correcting NR time step not to exceed final time, dt = " << scientific << newdt << ", gap = " << gap;
    }
	
    return newdt;
}






void newtonraphson :: info(std::ostream& of)
{
	stepsolver::info(of);
	
    of  << "\n Newton-Raphson non-linear solver"
	<< "\n Maximum number of iterations allowed    : " <<  maxIterations
	<< "\n Recompute tangent every # of iterations : " <<  recomputeTangent
	<< "\n Relative energy error for convergence   : " <<  ErelTolerance
	<< "\n Number of preiterations                 : " <<  preiterations
	<< "\n Maximum number of allowed trials        : " <<  maxNRtrials
	<< "\n Linesearch                              : " <<  (_linesearch ? "true" : "false");
}




void newtonraphson :: report(ostream& of)
{
	stepsolver::report(of);
	of  << "\n Total number of iterations     : " << totalIterations
		<< "\n Total number of solution steps : " << nsteps
		<< "\n Average iterations/step        : " << totalIterations/nsteps;
}


void newtonraphson :: solutionMessage(ostream &of, double mflops)
{
    double relatError = firstEerror == 0.0 ? 1.0 : Eerror/firstEerror;
    
    of  << "\n  Res: "  << scientific << setprecision(4) << Rerror 
		<< ", EErr: "	<< scientific << setprecision(4) << Eerror
		<< " ("			<< scientific << setprecision(1) << relatError 
		<< ")"
        << ", Ener: "   << scientific << showpos << setprecision(4) << theEnergy->energy
        << noshowpos;
	
    if (mflops > 0.0) of << ", Mflps: " << fixed << setw(6) << setprecision(2) << mflops;
    else              of << ", Mflps: --- ";
    of  << flush;
}






/* Solves one linear system of equations of the analysis K Du = R,
 where K is the global tangent, Du the increment in the solution
 and R the residual. First it assembles the residual and
 tangent.
 Returns 1 if ok, -1 if some problem
 */
bool newtonraphson :: solveIteration(std::ostream &of)
{
	DebugMessage("Solving NR iteration");
    double       tn, tn1, tna, alf;
	extern int global_iteration;
    
    // instant at which equilibrium is enforced
    tn1 = theIntegrator->theTime;
    tn  = tn1 - theTSControl->dt;
    alf = theIntegrator->getEqTime();
    tna = (1.0 - alf)*tn + alf*tn1;
    	
	
    // update tangent and residual
	if (iteration % recomputeTangent == 0 || totalIterations == 0)  
        _theAssembler->assembleResidualTangent(*theModel, tna);
	else														
        _theAssembler->assembleResidual(*theModel, tna);
	
	
    //  Solve the system of equations
	//  If using an iterative solver, there is no need to solve the system very accurately, in every iteration
	//  since eventually we will throw away the solution. We request only to reduce by 10000 the residual
	theLinearSolver->setRelTolerance(1e-4);
	bool ret = solveGlobalSystem(of);	
	theModel->incrementSolution(*theIntegrator);
	

//	theLinearSOE->getB().print();

	// update counters
	++totalIterations;
	++global_iteration;
	++iteration;
	++stepIters;
	
    return ret;
}





/* 
 Solve a Newton-Raphson loop, until convergence or error.
 return true if convergence ok, false if no convergence attained
 */
bool newtonraphson :: solveStep(std::ostream &of)
{
	DebugMessage("Solving NR step");
    extern int    global_iteration;
	extern double global_macheps;
	
    // we set the NR procedure as not converged initialize, and we increment the NR counter
    bool converged = false;
	iteration = 0;
    NRtrials++;
    
    // base the convergence on the relative energy error, relative to the first error
	// in some strange cases (only imposed bc with no loading and contact) the initial residual might be too low to
	// use as convergence criterion and it might be necessary to iterate a few times
	solveIteration(of);
	if (Eerror < 1e4*global_macheps)
		for (int i=0; i<preiterations; i++) solveIteration(of);
	
    firstResid  = Rerror;
    firstEerror = Eerror;
    double relatError  = 1.0;
	if (theTSControl->getName() == "Fixed stepsize") 	variableStepSize = false;
	else												variableStepSize = true;
	
    
	// first newton raphson loop trial.
	// the loop is ended because either: 1) it converges,   2) maxiter exceeded 
    bool   test(true);
    while (test)
    {
        solveIteration(of);
		
		relatError = firstEerror > 0.0 ? Eerror/firstEerror : 1.0;
        
        // convergence test. true if not convergence but can continue computing
        test =   relatError  >  ErelTolerance
			&&   Eerror      >  EabsTolerance
			&&   iteration   <  maxIterations
			&&   (isfinite(Eerror) != 0)
			&&   relatError	 <  divergenceRerror;
		
		// this serves to debug the login behing NR iterations, which is not as simple as one might think
		if (false)
		{
			of << endl;
			of << "True is    " << std::setw(9) << true << endl;
			of << "relatError " << setw(9) << relatError << "   ErelTolerance = " << ErelTolerance << ", test 1: " << (relatError  >  ErelTolerance) << endl;
			of << "Eerror     " << setw(9) << Eerror << ", test 2: " << (Eerror      >  ENERGY_ABS_TOLERANCE) << endl;
			of << "iteration  " << setw(9) << iteration << "  maxIterations " << maxIterations << ", test 3: " << (iteration   <  maxIterations) << endl;
			of << "finiteE    " << setw(9) << isfinite(Eerror) << "test 4: " << (isfinite(Eerror) != 0) << endl;
			of << "relatError " << setw(9) << relatError << " divergenceRerror " << divergenceRerror << ", test 5: " <<  (relatError <  divergenceRerror) << endl;
			of << "Test       " << setw(9) << test << endl;
		}
		
    }
    
    //	what to do with the newton raphson solution. 
	//  we evaluate it to see if it is suitable: 
	//     1) it has converged, 
	//     2) the estimated error is good enough
    
	
	// convergence because error is small
	if ( (relatError < ErelTolerance || Eerror < EabsTolerance) && !isnan(Eerror) )
	{
		// we only estimate the error in the integration when there is convergence
		// XXX EstimateIntegrationError(a->theIntegrator, a->theModel, theTSControl->dt, &(a->energies), a->theLinearSolver);
		converged = true;
	}
	
	
	// if fixed time stepping, and maxiter exceeded, there is nothing to be done
	else if ( !variableStepSize  && iteration == maxIterations)
	{
		of << "\n  <<<< WARNING: solution not converging. Newton-Raphson iterations exceeded." << endl;
		converged = false;
	}
	
	// if fixed time stepping, and the error is nan, there is nothing to be done
	else if ( !variableStepSize  && 
			 (iteration == maxIterations || relatError>=divergenceRerror || isnan(Eerror) || !isfinite(Eerror) )  )
	{
		of << "\n  <<<< WARNING: solution diverged. Fixed time stepping can not solve it." << endl;
		converged = false;
	}
	
	//  no convergence & maximum number of NR trials has been reached -> failure
    else if (variableStepSize  &&
			 (iteration    == maxIterations  || isnan(Eerror) || !isfinite(Eerror) || relatError > divergenceRerror) &&
             NRtrials      == maxNRtrials)
	{
		of << "\n  <<<< WARNING: solution not converging. Newton-Raphson trials exceeded." << endl;
        converged = false;
    }
	
	
	// there is no convergence and we have reached the minimum dt allowed -> failure
    else if (iteration == maxIterations && theTSControl->dt <= theTSControl->mindt)
	{
		of << "\n  <<<< WARNING: solution not converging. Time step size cannot be reduced further." << endl;
        converged = false;
	}
	
	
    // the analysis has not converged, but we can do something about it
	else
    {
        // halve the time step, but do not go below the minimum acceptable
        theTSControl->dt        = std::max<double>(0.5*theTSControl->dt, theTSControl->mindt);
        theIntegrator->theTime -= theTSControl->dt;
        iteration     = 0;                       
        global_iteration  = 0;
        of << "\n  <<<< WARNING: solution not converged. Reducing t-step size (trial #" << NRtrials+1 << "/" << maxNRtrials << ")" << endl;
        of << "\n  <<<< New time = " << scientific << setprecision(6) << theIntegrator->theTime 
		<< ", dt = " << theTSControl->dt << flush;
        
        // Go back to the start of the step, recovering the initial value of the variables
		theModel->rewindSolution();
        theIntegrator->startStep(*theModel, theTSControl->dt);
        
        // Then, reset the dU to zero, and introduce the imposed bc by modifying the dU of correspoding dofs
        theModel->resetSolutionIncrements();
		imposeConstrainedDofs();
        
        // Recursive call until problem solved or too many trials
        converged = solveStep(of);
    }
    
	
	/*  the newton raphson has converged, but we still need to check in some cases if the
	 (estimated) error in the solution of a transient problem is within the acceptable tolerance */
	if (theTSControl != NULL && typeid(*theTSControl) != typeid(fixedTSC) )
	{
		// if the time step has reached its minimum allowed value, do not even consider error estimation
		if (fabs(theTSControl->dt - theTSControl->mindt) <= global_macheps)
		{
			Message("\n  <<<< WARNING: integration error can not be further reduced. Mindt reached.");
			converged = false;
		}
		
		
		else
		{
			// compute the tolerance that the time step controller allows
			double tol = theTSControl->computeIntegrationTolerance(*theModel);
			theTSControl->setData(tol, iteration);
			
			// the error that the tsc has estimated is larger than the allowed tolerance, rewind and start a new step
			if (theTSControl->getError() >= 1.1*tol)
			{
				theIntegrator->theTime -= theTSControl->dt;
				theEnergy->increaseIntegrationError(- theEnergy->getStepError());
				theTSControl->dt = theTSControl->decreaseDt();
				
				iteration = global_iteration = NRtrials = 0;				
				theTSControl->setDt(theTSControl->dt);
				theIntegrator->theTime   += theTSControl->dt;
				
				of << "\n   <<<< WARNING: estimated error too large. Reducing dt.";
				of << "\n   New time = " << setprecision(5) << theIntegrator->theTime 
				<< ", dt = " << theTSControl->dt 
				<< flush;
				
				// Go back to the start of the step, recovering the initial value of the variables 
				theModel->rewindSolution();
				theIntegrator->startStep(*theModel, theTSControl->dt);
				
				// Then, reset the dU to zero, and introduce the imposed bc by modifying the dU of correspoding dofs
				theModel->resetSolutionIncrements();
				imposeConstrainedDofs();
				
				// Recursive call until problem solved or too many trials
				converged = solveStep(of);
			}
		} 
	}
	
    if (converged)
	{
		hmindt = std::min<double>(hmindt, theTSControl->dt);
        hmaxdt = std::max<double>(hmaxdt, theTSControl->dt);
	}
	
    return converged;
}
