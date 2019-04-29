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
 *  tscontrol.c
*
 */


#include "Analysis/TScontrol/tscontrol.h"

#include <iostream>
#include <iomanip>


#include <cstdlib>
#include <cstring>
#include <cmath>
#include <limits>
#include "Elements/eltype.h"
#include "Io/logger.h"
#include "Io/message.h"
#include "Io/usercommand.h"
#include "Math/linearsoe.h"
#include "Model/model.h"
#include "Analysis/Integrators/integrator.h"
#include "Analysis/energies.h"

#ifdef WITHMPI		
#include "mpi.h"
#endif


tscontrol :: tscontrol() :
	rejectedSteps(0),
	nSteps(1),
	userTolerance(0.01),
	tol(1.0), 	
	error(0.0),
	dt(0.0),
	mindt(-1.0),				// a negative number disables it
	maxdt(std::numeric_limits<double>::max()),  				// a large number disables it
	energy(0)
{
	errorHistory[0] = errorHistory[1] = 1.0;
	tolHistory[0]   = tolHistory[1]   = 1.0;
	dtHistory[0]    = dtHistory[1]    = 0.0;
}



tscontrol :: tscontrol(commandLine &cl) :
	rejectedSteps(0),
	nSteps(1),
	userTolerance(0.01),
	tol(1.0), 	
	error(0.0),
	dt(1.0),
	mindt(-1.0),				// a negative number disables it
	maxdt(std::numeric_limits<double>::max()),  				// a large number disables it
	energy(0)
{
	errorHistory[0] = errorHistory[1] = 1.0;
	tolHistory[0]   = tolHistory[1]   = 1.0;
	dtHistory[0]    = dtHistory[1]    = 0.0;
	
	for (int k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
		if       ( uc.keyword() == "max_dt" ) maxdt = uc.value();
		else if  ( uc.keyword() == "maxdt" )  maxdt = uc.value();
		else if  ( uc.keyword() == "min_dt" ) mindt = uc.value();
		else if  ( uc.keyword() == "dt" )     dt    = uc.value();
	}
}





double tscontrol :: adaptiveForwardDt()
{
	double newdt;
	
	newdt = chooseNewDt(dt, error, tol, dtHistory, errorHistory, tolHistory, parameter);
	
	
	// the time step should not be larger than the maximum nor smaller than the minimum
	if (newdt < mindt)
	{
		logger :: mainlog << "\n  Automatic time step " << setw(8) << scientific << setprecision(4) << newdt 
			<< " too small. Correcting to dt = " << mindt << flush;
		newdt = mindt;
	}
	else if (newdt > maxdt)
	{
		logger :: mainlog << "\n  Automatic time step " << setw(8) << scientific << setprecision(4)
			<< newdt << " too large. Correcting to dt = " << maxdt << flush;
		newdt = maxdt;
	}
	
	nSteps++;
    
  	return newdt;
}




void  tscontrol :: advanceDtVariables()
{
	dtHistory[1] = dtHistory[0];
	dtHistory[0] = dt;
}




void  tscontrol :: advance()
{
	double *tolP   = tolHistory;
	double *errorP = errorHistory;
	
	/* Copy tol, dt & error at t_{n+1} to t_n */
	tolP[1] = tolP[0];
	tolP[0] = tol;
	
	errorP[1] = errorP[0];
	errorP[0] = error;
}




double tscontrol :: computeIntegrationTolerance(model &m)
{
		
	double tolerance = 1e-16;
	
	if ( name == "Iteration based") return (double) tol;
	else
	{
		/* Compute the energy of the system loop over the elements 

		double systemEnergyNorm=0.0;
		element   el;
		energies elmtEnergies;
		int       r;
		
		el = m->firstelmt;
		while (el != NULL)
		{
			// get the energies in one element
			SetEnergiesToZero(&elmtEnergies);
			r = GetElementEnergies(el, &elmtEnergies);
			
			// and accumulate them into the global energies, in case something has been computed
			if (r>0) systemEnergyNorm += elmtEnergies.energy;
			el = el->next;
		}
		
		
		// Compute the tolerance
		tolerance = tsc->userTolerance*sqrt(systemEnergyNorm);
		*/
	}
	
	return tolerance;
}




double tscontrol :: decreaseDt()
{
	double ik, newdt;
	
	ik = 1.0/(orderMethod + 1.0);
	
	newdt = pow(tol/error, ik)*dt;
	rejectedSteps++;
	
	if (newdt < mindt) 	newdt = mindt;
	
 	return newdt;
}




bool tscontrol :: check(ostream& of)
{
	if (dt == 0.0) 
		of << "\n   Zero time step size in time stepsize control";
	else
		of << "\n   [ Stepsize control checked ]";
	return true;
}



/* Check if I, PI or PID time step controller has been selected and none error estimator has been chosen 
void      TSControlCheck(tscontrol tsc, dynestimator *dynest)
{
	// When tsc = NULL --> time-stepping is FIXED 
	if (tsc == NULL)	return;
	if (!IsIterationController(tsc) && dynest == NULL)
		ErrorMessage("\n An error estimator has to be selected to use with the %10s time step controller", tsc->name);
	
}
*/

	
void tscontrol :: info(ostream &of)
{
    of << "\n\n\n             T i m e    s t e p s i z e   s e l e c t i o n" << endl;
	of << endl << " Type : " << name;
}



/*  scan the command line
*  "stepping, type = [fixed|iteration| cfl], dt = xx, max_dt= xx, min_dt = xx, control = [i | pi | pid]"
*/
void tscontrol :: scan(commandLine &cl, tscontrol **t)
{
	tscontrol *tsc=0;
	
	for (int k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
		if       ( uc.keyword() == "type" && uc.option() == "fixed")	tsc = new fixedTSC(cl);
		else if  ( uc.keyword() == "type" && uc.option() == "automatic")tsc = new iterationTSC(cl);
		else if  ( uc.keyword() == "type" && uc.option() == "cfl")      tsc = new cflTSC(cl);
	}
	
	if (tsc != 0)
	{
		if (*t != 0) delete *t;
		*t = tsc;
	}
    
    if (*t == 0)
        *t = new fixedTSC(cl);
}




void tscontrol :: setLinks(energies &e, model &m)
{
	energy     = &e;
	linkedModel = &m;
}
	

void   tscontrol :: setData(double tol, int nriter)
{
	//energies->getStepError();
}



fixedTSC :: fixedTSC(commandLine &cl) :
	tscontrol(cl)
{
	name  = "Fixed stepsize";
	if (dt == 0.0) dt = 1.0;
}



fixedTSC :: fixedTSC(double dtc)
{
	name  = "Fixed stepsize with default parameters";
	dt = dtc;
}



double fixedTSC :: chooseNewDt(double dtc, double error, double tol, double *dtP, double *errorP, double *tolP, double *p)
{	
 	return dt;
}



void fixedTSC ::  info(ostream &of)
{
	tscontrol::info(of);
	of  << "\n No time step adaptation."
		<< "\n The first failure in a NR solution will cause the analysis to halt."
		<< "\n Time step size : " <<  setw(8) << scientific << dt << endl;
}



iterationTSC :: iterationTSC(commandLine &cl) :
	tscontrol(cl),
	targetIterations(6)
{
	name = "Iteration based";
	orderMethod = 1.0;
	if (dt == 0.0) dt = 1.0;
	
	for (int k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
		if       ( uc.keyword() == "target")	 targetIterations = static_cast<int>(uc.value());
	}
	
}



double iterationTSC :: chooseNewDt(double dt, double error, double tol, double *dtP, 
                                   double *errorP, double *tolP, double *p)
{
	double newdt;
	extern int global_iteration;
	
	if (global_iteration >= targetIterations*2)
		newdt = dt*0.5;
	else if(global_iteration >= targetIterations*1.2)
		newdt = dt*0.8;
	else if(global_iteration < targetIterations*0.8)
		newdt = dt*1.2;
	else if(global_iteration < targetIterations*0.5)
		newdt = dt*2.0;
	else
		newdt = dt;
	
 	return newdt;
}



void iterationTSC ::  info(ostream &of)
{
	tscontrol::info(of);
	of  << "\n\n Time step controller:"  << name
		<<   "\n Target number of iterations per step : " << targetIterations << endl;
}




iTSC :: iTSC(commandLine &cl) :
	tscontrol(cl)
{
	name = "Integral control";
	orderMethod = 2.0;
	if (dt == 0.0) dt = 1.0;
}


/* Strategy to select the new step siz:
dt_new = dt_old*(tol/error)^(1/p+1) ,
where 'p' is the order of the integration method  */
double iTSC :: chooseNewDt(double dt, double error, double tol, double *dtP, double *errorP, double *tolP, double *p)
{
	double newdt = pow(tol/error,p[0])*dt;
	return newdt;
}



void iTSC :: info(ostream &of)
{
	tscontrol::info(of);
	of  << "\n Control Parameter 1  :     " << p0 << endl;
}



cflTSC :: cflTSC(commandLine &cl) :
	tscontrol(cl), counter(0), computeEvery(1), scaletmin(1.0)
{
	name = "Time step selection based on CFL condition.";
	for (int k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
		if       ( uc.keyword() == "update" ) computeEvery = static_cast<int>(uc.value());
		else if  ( uc.keyword() == "scaledt") scaletmin  = uc.value();
	}
		
	if (dt == 0)
	{
		maxev = linkedModel->maxEigenvalueEstimate();
		dt    = 2.0/maxev;
	}
}



cflTSC :: cflTSC()
{
	name = "Time step selection based on CFL condition.";
}



/* Strategy to select the new step size:
*/
double cflTSC :: chooseNewDt(double dt, double error, double tol, double *dtP, double *errorP, double *tolP, double *p)
{
	double newdt;
	
	if ( counter++ % computeEvery == 0) 
	{
		maxev = linkedModel->maxEigenvalueEstimate() / scaletmin;
		if (maxev == 0.0) newdt = dtP[0];
		else 	newdt = 2.0/maxev;

#ifdef WITHMPI		
	// synchronize dt among computers
		double localdt = newdt;
		MPI::COMM_WORLD.Allreduce( &localdt, &newdt, 1, MPI::DOUBLE, MPI::MAX);
#endif
	}

	// repeat the same time step size
	else
		newdt = dtP[0];

	return newdt;
}



void cflTSC :: info(ostream &of)
{
	tscontrol::info(of);
	of  <<  "\n CFL condition for the second order differential equation."
	    <<  "\n Set for central differences integration."
	    <<	"\n Compute critical time step size every " << computeEvery
		<<  "\n Scale time step size by factor " << scaletmin << endl;
}



