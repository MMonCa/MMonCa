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
 *  stepsolver.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on Mon Dec 15 2003.
 *
 */

#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "Analysis/analysis.h"
#include "Analysis/assembler.h"
#include "Analysis/Dofs/dofset.h"
#include "Analysis/Loading/loading.h"

#include "Analysis/Stepsolvers/stepsolver.h"
#include "Analysis/Stepsolvers/linear.h"

#include "Main/feliks.h"
#include "Model/constraint.h"

#include "Io/logger.h"
#include "Io/message.h"
#include "Io/Postprocessor/postprocessor.h"
#include "Io/usercommand.h"

#ifdef WITHTBB
#include "tbb/task.h"
#include "tbb/task_scheduler_init.h"
#endif

#include <boost/foreach.hpp>



// a global variable indicating the counter for the predictor/corrector or multi-stage RK
int global_stage=0;

// a global variable that indicates if energy has to be recomputed or taken from the residual calculation
bool global_lazyenergy=false;


/* given a commandlist, find the command "stepsolver = xxx" and return the
 * corresponding stepsolver. If no solver is input return empty stepsolver
 */
stepsolver :: stepsolver(const commandLine& cl) :
    nsteps(0),
    stepIters(0),
	hmindt(HUGE),				// a large number initializes it
	hmaxdt(-1),					// a negative number initializes it
    cputime(0.0),
    NRtrials(0),
    maxNRtrials(MAX_NRTRIALS),
	variableStepSize(false),
	tmpVector(),
    _linesearch(false),
	theModel(0),
	theIntegrator(0),
	theLinearSOE(0),
	theLinearSolver(0),
	theTSControl(0),
	_theAssembler(0),
	theEnergy(0),
	stPostprocessor(0)
{
	for (size_t k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
		if		( uc.keyword() == "linesearch" && uc.option() == "on" )	_linesearch	= true;
	}			
}



/* given a commandlist, find the command "stepsolver = xxx" and return the
 * corresponding stepsolver. If no solver is input return empty stepsolver
 */
stepsolver :: stepsolver() :
	stepIters(0),
	nsteps(0),
	NRtrials(0),
	maxNRtrials(MAX_NRTRIALS),
	hmindt(HUGE),				// a large number initializes it
	hmaxdt(-1),					// a negative number initializes it
	variableStepSize(false),
	cputime(0.0),
	tmpVector(),
	theModel(0),
	theIntegrator(0),
	theLinearSOE(0),
	theLinearSolver(0),
	theTSControl(0),
	_theAssembler(0),
	theEnergy(0),
	stPostprocessor(0),
	_linesearch(false)
{
}


stepsolver :: ~stepsolver()
{
	if (_theAssembler != 0) delete _theAssembler;
	_theAssembler = 0;
}



#ifdef WITHTBB
// data structure to split advanceTime into several threads
class ioTask : public tbb::task{
	
public:
	ioTask(const stepsolver& theStepsolver, const model& theModel, 
           postprocessor* thePostprocessor, const energies& theEnergy, 
           const assembler& theAssembler, const double t) :
        _theStepsolver(theStepsolver), 
        _theMesh(theModel), 
        _thePostprocessor(thePostprocessor), 
        _theEnergy(theEnergy), 
        _theTime(t), 
        _theAssembler(theAssembler){}
	
	tbb::task* execute()
	{
		logger::logVariables(_theTime, _theMesh, _theEnergy, _theAssembler);		
		if (_thePostprocessor != NULL) _thePostprocessor->dumpSolution(_theMesh, _theTime);
		_theStepsolver.convergenceInfo();
		
		return NULL;
	}
	
private:
	const stepsolver&		_theStepsolver;
	const model&            _theMesh;
	postprocessor*			_thePostprocessor;
	const energies&         _theEnergy;
	const double			_theTime;
	const assembler&		_theAssembler;
};
#endif


/* Advances the global time counter for the integrator. This should be
 done after convergence is reached in the Newton-Raphson scheme. It also
 updates the values of the imposed bc.
 */
void stepsolver :: advanceTime(double finalTime, std::ostream &of)
{
	// when the step has finished, compute the energies
    if (theIntegrator->getEqTime() != 1.0) theModel->updateCurrentState(dofset::tn1);
    theModel->accumulateEnergy(*theEnergy, dofset::tn1);
	
	
#ifdef WITHTBBBB
	{
		// multithreaded version. It does not work. It executes the tasks sequentially
		tbb::task_scheduler_init init;
		
		// send io data
		ioTask   iTask(*this, *theModel, stPostprocessor, *theEnergy, theIntegrator->theTime, *_theAssembler);
		iTask.execute();
	}
#else
    // send information to file
	logger::logVariables(theIntegrator->theTime, *theModel, *theEnergy, *_theAssembler);
	if (stPostprocessor != 0) stPostprocessor->dumpSolution(*theModel, theIntegrator->theTime);
	convergenceInfo();
	
#endif
	
    // Compute new time step size
	extern double global_dt, global_dtn;
	global_dtn   = global_dt;
	double newdt = global_dt = adaptiveTimeStep(finalTime,of);
	theTSControl->setDt(newdt);
	
    
    bool ghostnds = false;
    if (theModel->recompute(ghostnds) == true) {
        
        
        if (theIntegrator->isImplicit())  // Implicit strategy of remeshing.
        {  
            if(ghostnds) //We have to mark the ghost nodes found and separate them from the problem
            {
                logger::mainlog << "\n\n WARNING: Ghost Nodes have appeared after remeshing!" << flush;
                logger:: sumf   << "\n\n WARNING: Ghost Nodes have appeared after remeshing!" << flush;
                cout            << "\n\n WARNING: Ghost Nodes have appeared after remeshing!" << endl;
            }
            
            //MMM We re-initialize the SOE, and after that we assign memory for the new tangent matrix with the reassignTangent function
            
            theLinearSolver->initialize(*theModel); // responsible for creating a linear SOE
            theLinearSOE = &(theLinearSolver->getTheLinearSOE()); 
            _theAssembler->reassignTangent(theLinearSOE->getA());      
            
        }
        
        else                              // Explicit strategy of remeshing.
            if(ghostnds) cout<< "Some ghostnodes have been assigned the gravity acceleration \n";
        
       /* messageAfterRemeshing(of);
        if (!solveStep(logger::mainlog) )
            {
                logger::mainlog << "\n\n WARNING:AFTER REMESHING ANALYSIS HAS NOT CONVERGED" << flush;
                logger:: sumf   << "\n\n WARNING:AFTER REMESHING ANALYSIS HAS NOT CONVERGED" << flush;
                cout            << "\n\n WARNING:AFTER REMESHING ANALYSIS HAS NOT CONVERGED" << endl;
            }*/
        
    }

    // Advance the degrees of freedom and rates according to integrator
    incrementStep();
    
    theIntegrator->advanceSolutionInTime(*theModel, theTSControl->dt);
	messageWhenAdvancingStep(of);
    
    // calculate predictors with the integrator and other factors
    theIntegrator->startStep(*theModel, theTSControl->dt);
	
    // update the imposed b.c. by changing the constrained dofs increments
    theModel->resetSolutionIncrements();
    imposeConstrainedDofs();
	
	// advance data in time step controller
	theTSControl->advance();
	
	// stepsolver info
	stepIters = 0;
}





bool stepsolver :: check(ostream &of)
{
	of  << "\n   [ Stepsolver checked ]";
	return true;
}



void stepsolver ::  convergenceInfo() const
{	
	if (theIntegrator->theTime > 0.0)
		logger::toConvFile(nsteps, theIntegrator->theTime, theTSControl->dt, stepIters, NRtrials, firstResid, Rerror);
}




/* Loop over all the constraints in the model, and set the nonzero
 degrees of freedom to their imposed value, as dictated by the proportionality 
 factor. 
 Depending on whether the algorithm is implicit or explicit, two different approaches
 are followed for imposing the boundary conditions:
 
 i) if the algorithm is explicit, the imposed dofs are directly stored in the
 nodal dofs.
 
 ii) if the algorithm is implicit, the imposed increments on the imposed b.c. are
 stored in the dof increment dU, and after the first iteration they will be
 added to the global dofs. This is crucial for plasticity/fracture problems.
 If the imposed dofs are applied directly to the affected nodes the deformation
 concentrates in adjacent elements and creates artificially large deformations
 on those areas, resulting in large gradients, hence plasticity, or damage,
 which will have to be reversed once the deformation is comunicated to the
 rest of the model ==> bad convergence of the NR solution.
 However, the present solution has a flaw: these "future" increments are
 stored in the variables delta, which really mean "the last increment".
 */
void  stepsolver :: imposeConstrainedDofs()
{
    // loop over the SP constraints in model
	list<constraint*>::iterator iter = theModel->spconstraints.begin();
    while (iter != theModel->spconstraints.end() )
    {
		(**iter).imposeConstrainedDOFs(theIntegrator->getTheTime(), theTSControl->dt, 
                                       theIntegrator->isImplicit(), *theIntegrator );
		++iter;
    }
    
    DebugMessage("Essential boundary conditions included.");
}




void stepsolver :: incrementStep()
{	
	extern int global_iteration;
	
	nsteps++;
	NRtrials     = 0;
	firstEerror  = 0.0;
    firstResid   = 0.0;
	global_iteration = 0;
}



// function to print the information of any integrator
void stepsolver :: info(ostream &of)
{
    of	<< "\n\n\n\n";
	printCentered(of, "S t e p   S o l v e r");
	of	<< "\n\n Solution strategy for each time step: " << name;	
}




// this is the function that is called when the information for the analysis is requested at the end of execution
void  stepsolver :: infoSummary(ostream &of)
{
	of <<  "\n\n The step solver:";
	of <<  "\n\tType                       : " <<  name;
	of <<  "\n\tTotal # of solution steps  : " << nsteps;
	of << endl;
}





bool stepsolver :: lineSearch(const longvector& initialResidual, longvector& theSolution)
{	
	// instant at which equilibrium is enforced
    double       tn, tn1, tna, alf;
    tn1 = theIntegrator->theTime;
    tn  = tn1 - theTSControl->dt;
    alf = theIntegrator->getEqTime();
    tna = (1.0 - alf)*tn + alf*tn1;	
	
	// The function we want to minimize through scaling is G(s) = delta * Residual[d_old + s*delta]
	// value of G for s=0.0;
	double G0 = theSolution.dot(initialResidual);
	linesearchMessage( 0.0, G0/G0, logger::mainlog);
	
	// store the (non-scaled) displacement in the temporary vector (just rename it to make more sense)
	longvector &DeltaD = tmpVector;
	DeltaD = theSolution;
	
	// value of G for s=1.0
	theModel->localizeSolutionIncrement(theSolution);
	theModel->incrementSolution(*theIntegrator);
	_theAssembler->assembleResidual(*theModel, tna);
	double G1 = DeltaD.dot( theLinearSOE->getB() );
	linesearchMessage(1.0, G1/G0, logger::mainlog);
	
	// undo the last increment in the solution
	theSolution = DeltaD; theSolution.changeSign();
	theModel->localizeSolutionIncrement(theSolution);
	theModel->incrementSolution(*theIntegrator);
	
	double stest=1.0;
	if ( G0*G1 < 0.0 )
	{
		// iterate by bisection trying to minimize G(s)
		double sleft = 0.0, sright = 1.0;
		double Gleft = G0 , Gright = G1, Gtest = G1;
		int lscounter = 0;
		while ( Gleft * Gright < 0.0 && fabs(Gtest/G0) > 0.2 && lscounter < 10)
		{
			// interpolation of (sl,Gl) --- (sr,Gr) by a line that should cross 0
			stest = sleft - Gleft*(sright-sleft)/(Gright-Gleft);
			
			// value of G for at the stest
			theSolution  = DeltaD; theSolution.scale( stest );
			theModel->localizeSolutionIncrement(theSolution);
			theModel->incrementSolution(*theIntegrator);
			_theAssembler->assembleResidual(*theModel, tna);
			Gtest = DeltaD.dot( theLinearSOE->getB() );
			linesearchMessage(stest, Gtest/G0, logger::mainlog);
			
			// undo the last increment in the solution
			theSolution = DeltaD;
			theSolution.scale( -stest );
			theModel->localizeSolutionIncrement(theSolution);
			theModel->incrementSolution(*theIntegrator);
			
			if (Gleft*Gtest > 0)
				sleft = stest, Gleft = Gtest;
			else
				sright = stest, Gright = Gtest;
			
			++lscounter;
		}
	}
	
	// set the final increments in the nodes
	theSolution  = DeltaD; 
	theSolution.scale( stest );
	
	return true;
}



void stepsolver :: linesearchMessage(const double s, const double Gs, std::ostream& of)
{
	of << endl << "       line search: s = " << scientific << setprecision(4) << s << ", G(s)/Go = " << setw(11) << scientific << setprecision(4) << Gs;
}



void stepsolver :: messageWhenAdvancingStep(std::ostream &of)
{
	of  << "\n\nSolution at (pseudo) time = " << scientific << setprecision(6) << theIntegrator->theTime
		<< ", dt = "                      << theTSControl->dt
		<< ", Step = "                    << fixed << setprecision(0) << nsteps;	
	factorcombo::printAllValues(theIntegrator->theTime, of);
}



void stepsolver :: messageAfterRemeshing(std::ostream &of)
{
	of  << "\n\n Solution without time-advance after remeshing: "
    <<"\nSolution at (pseudo) time = " << scientific << setprecision(6) << theIntegrator->theTime
    << ", dt = "                      << theTSControl->dt
    << ", Step = "                    << fixed << setprecision(0) << nsteps;	
	factorcombo::printAllValues(theIntegrator->theTime, of);
}




void stepsolver :: report(std::ostream &of)
{
	of  << "\n\n"
        << "\n----------------------------------------------------------------"
        << "\n\t\t\t Solver Report" 
        << "\n----------------------------------------------------------------"
        << "\n\n Total number of solution steps : " << nsteps
		<<   "\n Total time in solver           : " << cputime;
}




// the variable ss contains a step solver, defined by default in case no other is provided.
// if the cl defines a new stepsolver, the one previously in ss must be deallocated and
// replaced
void stepsolver :: scan(commandLine &cl, stepsolver **ss)
{
	stepsolver *ns=0;
	for (int k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
        if       ( uc.keyword() == "type" && uc.option() == "linear")	 ns = new linearss(cl);
        else if  ( uc.keyword() == "type" && uc.option() == "newton")    ns = new newtonraphson(cl);
	}	
	
	if (ns != 0)
	{
		if (*ss != 0) delete *ss;
		*ss = ns;
	}
}





void stepsolver :: setLinks(model &m, integrator &i, linearSOE &soe, linearSolver &ls, 
							tscontrol &tsc, energies &e, postprocessor &pp)
{
	theModel        = &m;
	theIntegrator   = &i;
	theLinearSOE    = &soe;
	theLinearSolver = &ls;
	theTSControl    = &tsc;
	theEnergy       = &e;
	stPostprocessor = &pp;
	
#ifdef WITHMPI 
	if(i.isImplicit() && _theAssembler == 0)
		_theAssembler   = new MPIImplicitAssembler( soe.getB(), soe.getA() );
	else	
		dynamic_cast<MPIExplicitAssembler*>(_theAssembler)->setLinkedMesh(m);

#else
	if (i.isImplicit() && _theAssembler == 0)
		_theAssembler   = new assembler( soe.getB(), soe.getA() );
#endif
	
}




/* Solve the global finite element system of equations, including line search
 */
bool stepsolver :: solveGlobalSystem(std::ostream &of)
{
	// instant at which equilibrium is enforced
    double       tn, tn1, alf;
    tn1 = theIntegrator->theTime;
    tn  = tn1 - theTSControl->dt;
    alf = theIntegrator->getEqTime();
	
	
    // this is a temporary vector that we use to compute the energy residual norm.
	// To avoid having to allocate memory for it every time this routine is called
	// we do it once and check to see if it is big enough every time
	
	// store the initial residual
	longvector& initialResidual = tmpVector;
	initialResidual = theLinearSOE->getB();
	
	
    // Solve the global system K dU = R, where K may have dynamic contributions,
	// as defined by each integrator. After the solution, theResidual is
	// overwritten by the dof increment dU. 
    double   cpu1, mflops;

//    theLinearSOE->getA().print();
//    theLinearSOE->getB().print();
    bool ret  = theLinearSolver->solve(*theLinearSOE, cpu1, mflops);
//	theLinearSOE->getB().print();

    cputime  += cpu1;
	
    // the error in the energy norm is obtained as the dot product of the residual
	// (before the solution of the equations) and the "displacement" increment
	// as a result of this residual
    longvector &theSolution = theLinearSOE->getB();
    Eerror = fabs( theSolution.dot( initialResidual) );
    Rerror = initialResidual.norm();
    theModel->accumulateEnergy(*theEnergy, dofset::tn1);
    solutionMessage(of, mflops);
	
	if (_linesearch) lineSearch(initialResidual, theSolution);
	theModel->localizeSolutionIncrement(theSolution);
	theSolution.setZero();
	
    return ret;
}




/* Solve the global finite element system of equations, which must be already
 built. After the solution phase, the nodal dofs must be updated and the
 residual cleared.
 return 1 if ok, -1 if the system can not be solved
 
 bool stepsolver :: solveGlobalSystem(std::ostream &of)
 {
 // this is a temporary vector that we use to compute the energy residual norm.
 // To avoid having to allocate memory for it every time this routine is called
 // we do it once and check to see if it is big enough every time
 tmpVector = theLinearSOE->getB();
 
 // Solve the global system K dU = R, where K may have dynamic contributions,
 //    as defined by each integrator. After the solution, theResidual is
 //    overwritten by the dof increment dU. If explicit, no K is needed
 double   cputime, mflops;
 
 bool ret  = theLinearSolver->solve(*theLinearSOE, cputime, mflops);
 Mflops    = (Mflops * total_iter + mflops) / (total_iter + 1);
 cputime  += cputime;
 
 
 // the error in the energy norm is obtained as the dot product of the residual
 //    (before the solution of the equations) and the "displacement" increment
 //    as a result of this residual
 longvector &theSolution = theLinearSOE->getB();
 Eerror = fabs( theSolution.dot(tmpVector) );
 Rerror = tmpVector.norm();
 solutionMessage(of, mflops);
 
 // the solution from the linear system of eqs in contained in theResidual, and
 // it is transferred next to the nodes
 theModel->localizeSolutionIncrement(theSolution);
 theSolution.setZero();
 
 // update counters
 incrementIterations();
 
 return ret;
 }
 */




bool stepsolver ::	solveIteration(std::ostream &of)
{
	of << std::endl << "This type of stepsolver does not need to iterate." << std::flush;
	return true;
}

