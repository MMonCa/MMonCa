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
 *  linearss.cpp
 *  MacFELIKS
 *
 *  Created by Ignacio Romero on Thu Sep 11 2003.
 *  Copyright (c) 2003 Universidad Politecnica de Madrid. All rights reserved.
 *
 */

#include "linear.h"

#include <cmath>
#include <iostream>
#include <iomanip>


#include "Analysis/analysis.h"
#include "Model/model.h"
#include "Io/logger.h"
#include "Analysis/assembler.h"
#include "Io/message.h"
#include "Io/message.h"
#include "Main/feliks.h"
#include "Math/linearsoe.h"
#include "Math/Sparse/sparse.h"




linearss :: linearss() :
	Etolerance(ENERGY_TOLERANCE),
	Rtolerance(RESIDUAL_TOLERANCE)
{
	implicit    = true;
	nsteps      = 0;
	name        = "Linear";
}


linearss :: linearss(commandLine &cl) :
	Etolerance(ENERGY_TOLERANCE),
	Rtolerance(RESIDUAL_TOLERANCE)
{
	implicit    = true;
	nsteps      = 0;
	name        = "Linear";
}



linearss :: ~linearss()
{
}




double linearss :: adaptiveTimeStep(double finalTime, std::ostream &of)
{
	double newdt, t = theIntegrator->theTime, dt = theTSControl->dt;
	extern double global_macheps;
	
	theTSControl->advanceDtVariables();

	if (t == 0.0 && dt == 0.0)
	{
		newdt = 1000.0/theModel->maxEigenvalueEstimate(); // the 1000 factor is a guess
		of << "\n  No starting time step size. Guessing with "<< newdt << flush;
	}
    else 		
		newdt = theTSControl->adaptiveForwardDt();
    
	
	// final check. Do not go over final time
	// some extra space to allow for fp errors
	double gap = (finalTime-t) + finalTime*0.5*(nsteps+1)*global_macheps;
	if (finalTime != -1 && newdt-gap > 0.0)
    {
		of << endl;
		of << "\n  Correcting time step not to exceed final time, dt = " << newdt << ", gap = " << gap;
        newdt = finalTime - theIntegrator->theTime;
		of << "\n  Correcting time step not to exceed final time, dt = " << newdt << ", gap = " << gap;
    }
		
    return newdt;
}




void linearss :: info(ostream &of)
{
	stepsolver::info(of);
    of  << "\n Step solver for linear problems.";
	of  << "\n Only one linear system of equations is solved per time step"
		<< "\n The tangent is reused if possible" << endl;
}




void linearss :: solutionMessage(ostream &of, double mflops)
{
    double relatError = firstEerror == 0.0 ? 1.0 : Eerror/firstEerror;
    
    of  << endl 
		<< "  Res: "  << scientific << setprecision(6) << Rerror 
		<< ", EErr: " << scientific << setw(8) << Eerror
		<< " ("       << scientific << setprecision(1) << relatError 
		<< ")";
	
    // approximate condition number
	/*
    m     = GetAInLinearSOE(a->theLinearSOE);
    cond  = SparseMatrixCondition(m);
    if (cond   != -1) Message("%2.1e", cond);
    else            Message("     ");
    */
	
    if (mflops > 0.0) of << ", Mflps: " << fixed << setprecision(2) << mflops;
    else              of << ", Mflps: --- ";
}



/* Solves one linear system of equations of the analysis K Du = R,
where K is the global tangent, Du the increment in the solution
and R the residual. First it assembles the residual and
tangent.
Returns 1 if ok, -1 if some problem
*/
bool linearss :: solveStep(std::ostream &of)
{
    // instant at which equilibrium is enforced
    double tn1 = theIntegrator->theTime;
    double tn  = tn1 - theTSControl->dt;
    double alf = theIntegrator->getEqTime();
    double tna = (1.0 - alf)*tn + alf*tn1;
    
    // update tangent and residual only if time step size changes
	if (nsteps == 1 || theTSControl->dt != theTSControl->getPreviousDt() )
		_theAssembler->assembleResidualTangent(*theModel, tna);
	else
		//_theAssembler->assembleResidual(*theModel, tna, 'A');
		_theAssembler->assembleResidualTangent(*theModel, tna);
	
    //  Solve the system of equations
	//  If using an iterative solver, there is no need to solve the system very accurately, in every iteration
	//  since eventually we will through away the solution. We request only to reduce by 1000 the residual
	theLinearSolver->setRelTolerance(1e-12);
    bool ret = solveGlobalSystem(of);
	theModel->incrementSolution(*theIntegrator);
	
	// update counters
	extern int global_iteration;
	++global_iteration;
	++stepIters;
	
	
	// a check. Uncomment to see that the tangent is correct and the residual after a solve is zero
	// normally, leave this commented. The residual should be zero, and it is a waste of CPU time to compute it
	_theAssembler->assembleResidual(*theModel, tna, 'A');
	longvector& res = theLinearSOE->getB();
	of << "\n  Res: " << scientific << setprecision(6) << res.norm() << flush;
	
    return ret;
}


