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
 *  festatic.c
 *  feliks
 *
 *  Created by Ignacio Romero on 20/8/05.
 *
 *
 */

#include "Analysis/assembler.h"
#include "Analysis/Dofs/dofset.h"
#include "Analysis/TScontrol/tscontrol.h"
#include "Analysis/FEanalysis/festatic.h"
#include "Analysis/Integrators/integrator.h"
#include "Analysis/Integrators/quasistatic.h"
#include "Analysis/Stepsolvers/newtonraphson.h"
#include "Analysis/Stepsolvers/linear.h"
#include "Analysis/Loading/loading.h"
#include "Analysis/Loading/surfaceloading.h"

#include "Main/feliksinterface.h"

#include "Io/io.h"

#include "Math/LinearSolvers/direct.h"
#include <cstring>
#include <boost/foreach.hpp>

using namespace std;

extern ofstream       mainlog;





staticFEanalysis :: staticFEanalysis(const feliks::mmoncaInterface& mi) :
	theIntegrator(0),
    finalTime(1.0),
    nsteps(0)
{
	theIntegrator   = new quasiStatic();
	theStepsolver   = new linearss();
	theLinearSolver = new direct();
	theTSControl    = new fixedTSC();	
	
	createFEdata(mi);
	theModel->fill(mi);
	theModel->initialize();
    dofs_per_node = theModel->getDofsPerNode();
    
	theLinearSolver->initialize(*theModel); // responsible for creating a linear SOE

	theLinearSOE = &(theLinearSolver->getTheLinearSOE());
	theStepsolver->setLinks(*theModel,
                            *theIntegrator, *theLinearSOE, *theLinearSolver, *theTSControl, theEnergies, *thePost);

	theTSControl->setLinks(theEnergies, *theModel);
	FEanalysis::check(logger::mainlog);

	info(logger::mainlog);
	theModel->info(logger::mainlog);
	theModel->print("all", logger::mainlog);
	infoSummary(logger::sumf);
}






staticFEanalysis :: staticFEanalysis() :
	finalTime(1.0),
	nsteps(0)
{
	// default integrator
	theIntegrator = new quasiStatic();
	
	extern analysisTypeT global_analysistype;
	global_analysistype = FEANALYSIS_STATIC;
}


staticFEanalysis :: ~staticFEanalysis()
{
	delete theIntegrator;
}



bool staticFEanalysis :: check()
{
	bool ret = true;

	FEanalysis :: check();
	theIntegrator->check();
        
	return ret;
}




void staticFEanalysis :: info(ostream &of)
{
	FEanalysis :: info(of);
	of << "\n   Final time          : " << finalTime;

    theIntegrator->info(of);
}





bool staticFEanalysis :: specificSolve()
{
	bool ret(true);
	
	// incremental solution
    bool converged=true;
    extern double global_macheps;
    
	// initialize static solution
	(theStepsolver->getAssembler()).setNumberOfTimeDerivatives(0);
	theIntegrator->initialize(*theModel, *theLinearSolver, *theLinearSOE, theStepsolver->getAssembler() );
	
	
    // advance until 'finaltime'
	while (finalTime - theIntegrator->theTime > (0.5*finalTime* nsteps*global_macheps) && converged == true)
	{
		// break the solution loop if the there is an error in the NR solution
		if ( !converged ) break;

		theStepsolver->advanceTime(finalTime, logger::mainlog);
		converged = theStepsolver->solveStep(logger::mainlog);
		++nsteps;
	}
	
	if ( !converged )
	{
		logger::mainlog << "\n\n WARNING: ANALYSIS HAS NOT CONVERGED" << endl;
		logger::sumf    << "\n\n WARNING: ANALYSIS HAS NOT CONVERGED" << endl;
		cout            << "\n\n WARNING: ANALYSIS HAS NOT CONVERGED" << endl;
		ret = false;
	}

	// before quitting, dump the data to the files
	theStepsolver->convergenceInfo();
    
    if (theIntegrator->getEqTime() != 1.0) theModel->updateCurrentState(dofset::tn1);
    theModel->accumulateEnergy(theEnergies, dofset::tn1);
	logger::logVariables(theIntegrator->getTheTime(), *theModel, theEnergies, theStepsolver->getAssembler() );
	FEanalysis :: info(logger::sumf);
	if (thePost.get() != NULL) thePost->dumpSolution(*theModel, theIntegrator->theTime);

	FEanalysis::finalMessage(logger::mainlog, ret);
	return ret;
}
