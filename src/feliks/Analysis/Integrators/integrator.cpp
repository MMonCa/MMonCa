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
 * integrator.cpp
 *
 * ignacio romero
 * june 2000, converted to C++ in dec 2005
 *
 */


#include <cstdlib>
#include <cstring>
#include <cmath>
#include <pthread.h>

#include "Analysis/assembler.h"
#include "Main/feliks.h"
#include "Elements/eltype.h"
#include "Io/logger.h"
#include "Io/usercommand.h"
#include "Math/linearsoe.h"
#include "Model/model.h"

#include "Analysis/Integrators/quasistatic.h"


#ifdef WITHTBB
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#endif


/* this is an ugly thing we are forced to use, for the moment: a global variable with
information about the integration scheme. He have to do this because there is certain
data from the integration that needs to be passed to the element, for the tangent
computation. This should not be like this, since the tangent should be integratior
independent. It is like that for most cases, but in the case of LMS integrators a big
saving is obtained by noting that the velocity and acceleration increments are
proportional to the displacement increments, the main variable.
In a similar way, for a explicit computation the inertial terms should not be included
in the residual, since they are actually part of the unknowns.
*/
double ftgstiff;	// this factor multiplies the stiffness matrix
double ftgdamp; 	// this factor multiplies the damping matrix ... in 2 order eq.
double ftgmass;		// this factor multiplies the mass matrix ... in 2 order eq






integrator :: integrator() :
	theTime(0.0),
	integrationOrder(0),
	nodallyBased(true),
	variableMass(false),
	spatialUpdates(false)
{
}


integrator :: integrator(const commandLine& cl) :
	theTime(0.0),
	integrationOrder(0),
	nodallyBased(true),
	variableMass(false),
	spatialUpdates(false)
{
	for (int k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
		if		( uc.keyword() == "variablemass" ) variableMass = true;
	}	
	
}




integrator :: ~integrator()
{	
}




// These are generic functions to advance dofsets. Specific integrators can define
// their own updates to gain efficiency (if not all the rates are used) or to
// do more sophisticated operations.

void integrator :: advanceDofset(vectorDofset &UDS) const
{
	UDS.displacement(dofset::tn) = UDS.displacement() = UDS.displacement(dofset::tn1);
	UDS.velocity(dofset::tn)     = UDS.velocity()     = UDS.velocity(dofset::tn1);
 	UDS.acceleration(dofset::tn) = UDS.acceleration() = UDS.acceleration(dofset::tn1);		
}


void integrator :: advanceDofset(scalarDofset &PDS) const
{
	PDS.displacement(dofset::tn) = PDS.displacement() = PDS.displacement(dofset::tn1);
	PDS.velocity(dofset::tn)     = PDS.velocity()     = PDS.velocity(dofset::tn1);
 	PDS.acceleration(dofset::tn) = PDS.acceleration() = PDS.acceleration(dofset::tn1);
}



void integrator :: advanceDofset(rotationDofset &DS) const
{    
	DS.rotation(dofset::tn)				= DS.rotation()            = DS.rotation(dofset::tn1);
	DS.spatialVelocity(dofset::tn)		= DS.spatialVelocity()     = DS.spatialVelocity(dofset::tn1);
	DS.spatialAcceleration(dofset::tn)  = DS.spatialAcceleration() = DS.spatialAcceleration(dofset::tn1);
	DS.spatialStepRotationVector().setZero();
}



void integrator :: advanceDofset(directorDofset &DS) const
{    
	DS.bodyVelocity(dofset::tn)     = DS.bodyVelocity()     = DS.bodyVelocity(dofset::tn1);
	DS.bodyAcceleration(dofset::tn) = DS.bodyAcceleration() = DS.bodyAcceleration(dofset::tn1);
	DS.rotation(dofset::tn)         = DS.rotation()         = DS.rotation(dofset::tn1);
	DS.bodyStepRotationVector().setZero();	
}





// we call this function when a solution has been reached that we take as
// good. All the counters are advanced, and a message is passed to all
// the components of the analysis that the solution is finally good,
// and they can update their data
void integrator :: advanceSolutionInTime(model &m, const double dt)
{
	advanceTimeCounters(dt);
	m.advanceInTime(*this, dt);
}



void integrator :: advanceTimeCounters(const double dt)
{
	extern double global_tn, global_tna, global_tn1, global_dt, global_dtn;
    
    global_tn   = theTime;
	theTime    += dt;
    global_tn1  = theTime;
    global_tna  = (1.0 - eqtime)*global_tn + eqtime*global_tn1;
	global_dtn  = global_dt;
    global_dt   = dt;	
}




bool integrator :: check()
{
	bool ret = true;	
	return ret;
}




// most of the integrators need not implement this function
void integrator :: incrementElementSolution(element& e) const
{
}



#ifdef WITHTBB
// data structure for tbb threaded version
class tbbIncrementSolution
{
public:
	tbbIncrementSolution(const integrator& integr, vector<node*>& nds) : 
		my_nodes(&nds),
		my_integrator(&integr) {}
	
	void operator()( const tbb::blocked_range<size_t>& r) const
	{
		vector<node*>       &theNodes( *my_nodes );
		const integrator    &theIntegrator(*my_integrator);
		
		for (size_t a=r.begin(); a!=r.end(); a++)	
			theNodes[a]->incrementSolution(theIntegrator);
	}
	
private:
	vector<node*>		*const my_nodes;
	const integrator    *const my_integrator;
};
#endif


/* from the solution phase, the nodes contain in dU the increment for each
degree of freedom. In this function we accumulate these values in
the nodal dofs and clear the dU, o
*/
void integrator :: incrementSolution(vector<node*> &nodes) const
{
#ifdef WITHTBB
	// tbb threaded version with fixed n of threads, and with automatic number
	tbb::parallel_for(tbb::blocked_range<size_t>(0,nodes.size()), tbbIncrementSolution(*this, nodes) );	
	
#else
	//serial version
	vector<node*>::iterator iter = nodes.begin();
	while (iter != nodes.end()  )
	{
		(*iter)->incrementSolution(*this);
		++iter;
	}
#endif
}


/* from the solution phase, the nodes contain in dU the increment for each
degree of freedom. In this function we accumulate these values in
the nodal dofs and clear the dU, o
*/
void integrator :: incrementSolution(vector<element*> &elements) const
{
	vector<element*>::iterator iter = elements.begin();
	while (iter != elements.end()  )
	{
		incrementElementSolution( **iter );
		++iter;
	}
}




void integrator :: initialize(model &m, linearSolver &linsolver, linearSOE &theSOE, assembler& as)
{		
	// set up initial conditions and rates.
    m.imposeInitialConditions();
    m.imposeInitialRates();		
}




// function to print the information of any integrator
void integrator :: info(ostream &of)
{
	of << "\n\n";
	printCentered(of, "I n t e g r a t o r");
	of << "\n";
}



// Replace whatever comes in *i by the integrator defined in commandline, if a
// correct integrator is defined. Otherwise, leave as it is
void integrator :: scan(commandLine &cl, integrator **i)
{
	integrator *ni=0;
	for (int k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
		if       ( uc.keyword() == "type" && uc.option() == "quasistatic")	ni = new quasiStatic(cl);
	}
	
	if (ni != 0)
	{
		if (*i != 0) delete *i;
		*i = ni;
	}
}

