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
 * analysis.cpp
 *
 *  feliks
 *
 *  Created by Ignacio Romero on Nov 11 1999.
 *  Copyright (c) 2003 Universidad Politecnica de Madrid. All rights reserved.
 *  converted to c++ by i.r. in march 2006
 * 
 * defines properties of general type of finite elmement analyses. 
 */

#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <typeinfo>

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <memory>
#include <unistd.h>


#include "Main/feliksinterface.h"
#include "Analysis/analysis.h"
#include "Analysis/assembler.h"
#include "Analysis/Propfactors/propfactor.h"


// FEanalysis
#include "Analysis/FEanalysis/festatic.h"
#include "Analysis/Integrators/quasistatic.h"

#include "Analysis/Loading/pointload.h"
#include "Analysis/Loading/surfaceloading.h"
#include "Analysis/Loading/nodesetloading.h"
#include "Analysis/Loading/volumeloading.h"


//stepsolvers
#include "Analysis/Stepsolvers/newtonraphson.h"

//linearsolvers
#include "Math/LinearSolvers/linearsolver.h"
#include "Math/LinearSolvers/direct.h"

//time step controler
#include "TScontrol/tscontrol.h"

// model related
#include "Elements/element.h"
#include "Model/constraint.h"
#include "Model/initrate.h"
#include "Model/model.h"

#include "Materials/material.h"
#include "Materials/Smallstrain/elastic.h"
#include "Materials/Smallstrain/alloy.h"


#include "Elements/eltype.h"
#include "Elements/Solid/Smallstrain/solid.h"
#include "Io/message.h"
#include "Io/logger.h"
#include "Io/usercommand.h"
#include "Main/feliks.h"
#include "Math/linearsoe.h"
#include "Math/vector.h"
#include "Math/tensor.h"
#include "Math/Sparse/sparse.h"

// postprocessor
#include "Io/Postprocessor/postprocessor.h"


#define DEBUG_SCANNER 0

using namespace std;

extern ofstream       mainlog;


/*  global variables that might be used elsewhere
    I am aware global variables are generally discouraged. 
	However, some information about the overall analysis needs to be obtained,
	from time to time, in different parts of the analysis. 	It would be 
	inconvenient to pass the whole analysis to every function, so part of 
	this global information is put into global variables. 
*/
double global_tn, global_tna, global_tn1, global_dt, global_dtn;
double global_error;
double global_wavenumber;
int    global_iteration;                    // iteration within each NR loop
int    global_complexdata=0;				// complex degrees of freedom
analysisTypeT    global_analysistype;



FEanalysis :: FEanalysis() :
	theModel(), 
	theLinearSOE(0),
    theStepsolver(0),
    theLinearSolver(0),
	theTSControl(0),
	mass(0),
    thePost(0),
	dofs_per_node(0),
	interactive(false)
{
	theModel = new model();
	theEnergies.initialize();		
}



FEanalysis :: ~FEanalysis()
{
	if (theLinearSolver!=0) delete theLinearSolver;
	if (theTSControl  != 0) delete theTSControl;
	if (theStepsolver != 0) delete theStepsolver;
	if (theModel      != 0) delete theModel;
}



/* a set of checks to be done before starting the analysis. This is not
complete yet. More checks need to be added
*/
bool FEanalysis :: check(ostream &of)
{    
	of << "\n\n\n Checking analysis and model ...";
    bool ok = theModel->check(of);

	if (theStepsolver != 0) ok = ok && theStepsolver->check(of);
	if (theTSControl  != 0) ok = ok && theTSControl->check(of);
	
	logger::checkLoggers(*theModel);
    ok = ok && specificCheck(of);
	if (!ok) ErrorMessage("Error in the analysis. Quitting FELIKS");

    return ok;
}




// now, with all the data in the analysisDescription, create FE structure
// in the required order, independently of the order provided by the user
void FEanalysis :: createFEdata(const feliks::mmoncaInterface& mi)
{
    
    // create pure materials
    std::vector<feliks::mmoncaMaterial>::const_iterator iter = mi.theMaterials_.begin();
    std::vector<elasticMaterial*> thePureMaterials;
    size_t a=0;
    while ( iter != mi.theMaterials_.end() )
    {
        material* mat = new elasticMaterial((*iter).name_, (*iter).ref_temp_, (*iter).thexpansion_,
                                            (*iter).E_, (*iter).nu_, (*iter).eigenstrains_);
        material::allMaterials[a++] = mat;
        thePureMaterials.push_back(dynamic_cast<elasticMaterial*>(mat));
        ++iter;
    }
    
    // create alloy
    smallStrainMaterial* alloy = new alloyMaterial("the alloy", thePureMaterials);

    
    // create the unique element type, which must be hold the alloy
    eltype*  se  = new solid(*alloy);
    eltype::allEltypes[0] = se;
}




void FEanalysis :: finalMessage(ostream &of, bool ok)
{
	of  << "\n\n\n End of computations.";
	if (ok)  of << "\n Analysis finished correctly." << endl;
	else     of << "\n The analysis could not be completed." << endl;
	
	if (theStepsolver != 0) theStepsolver->report(of);
}




void FEanalysis :: info(ostream &of)
{
	time_t walltime = time(NULL);

    of << "\n\n\n\n";
	printCentered(of, "A n a l y s i s    i n f o r m a t i o n");	
	of << "\n\n General properties:";
	of << "\n   Started on          : " << ctime(&walltime);
	of <<   "   Username            : " << getlogin();

    char host[256];
    if (gethostname( host,  256) == 0)
        of << "\n   Hostname            : " << host;


	if (theStepsolver != 0) theStepsolver->info(of);
	if (theTSControl  != 0) theTSControl->info(of);
	if (theLinearSolver!=0) theLinearSolver->info(of);
	eltype         :: listEltypes(of);
	material       :: listMaterials(of);
	factorcombo    :: listPropfactorCombinations(of);
	logger         :: listLoggers(of);
	if (thePost.get() != 0) thePost->print(of);	
}




void FEanalysis :: infoSummary(ostream &of)
{
	of << "\n\n                A n a l y s i s    i n f o r m a t i o n";
	
	// general information
	of << "\n\n General properties:";
	
	
	//if (theIntegrator != 0) theIntegrator->info(of);
	if (theStepsolver != 0) theStepsolver->infoSummary(of);
	if (theLinearSolver!=0) theLinearSolver->info(of);
	if (theLinearSOE  != 0) theLinearSOE->info(of);	
	if (theTSControl  != 0) theTSControl->info(of);
}





void FEanalysis :: setInteractive()
{
	interactive = true;
}




/* this is the main solution function, the one that drives the solution of
the problem. If the solution is interactive though, the flow is changed
to the user.
*/
bool  FEanalysis :: solve()
{
	logger::startAnalysisTimer();
	bool ret(true);
	
	logger::mainlog
		<< "\n\n\n\n"
		<< "\n\t            S t a r t i n g    A n a l y s i s"
		<< "\n\t            ----------------------------------"
		<< "\n";
	
	
	initialize();
    ret = specificSolve();

	return ret;
} 

