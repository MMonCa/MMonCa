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
 * quasistatic.cpp
 *
 * i. romero july 2000
 *
 * quasistatic integrator
 *
 * the velocity is stored as the increment of the displacements over the pseudotime. This false velocity
 * helps in contact detections
 */

#include <cstdio>
#include <cstring>

#include "Analysis/Integrators/quasistatic.h"
#include "Analysis/Integrators/integrator.h"
#include "Analysis/assembler.h"
#include "Elements/eltype.h"
#include "Io/message.h"
#include "Io/usercommand.h"
#include "Math/linearsoe.h"
#include "Math/tensor.h"
#include "Math/LinearSolvers/linearsolver.h"
#include "Model/model.h"
#include "Model/Node/node.h"


using namespace blue;

quasiStatic :: quasiStatic()
{
	name = "Quasistatic";
	implicit = true;
	eqtime   = 1.0;
	integrationOrder = 0;
	
	spatialUpdates = true;
}
	


quasiStatic :: quasiStatic(const commandLine &cl)
{
	name = "Quasistatic";
	implicit = true;
	eqtime   = 1.0;
	integrationOrder = 0;
	spatialUpdates = true;
	
	
	for (size_t k=1; k< cl.size(); k++)
	{
		const usercommand  &uc = cl[k];
		if  ( uc.keyword() == "rotupdates" && uc.option() == "spatial" ) spatialUpdates = true;
	}		
}



quasiStatic :: ~quasiStatic()
{
}



/* for the quasi-static integrator there is only one thing to be done, copy
the configuration variables at time tn+1 to the slot at tn
*/
void quasiStatic :: advanceDofset(vectorDofset &DS) const
{    
	DS.displacement(dofset::tn) = DS.displacement() = DS.displacement(dofset::tn1);
}


void quasiStatic :: advanceDofset(scalarDofset &DS) const
{    
	DS.displacement(dofset::tn) = DS.displacement() = DS.displacement(dofset::tn1);
}


void quasiStatic :: advanceDofset(rotationDofset &DS) const
{    
	DS.rotation(dofset::tn)	= DS.rotation() = DS.rotation(dofset::tn1);
	DS.spatialStepRotationVector().setZero();
}


void quasiStatic :: advanceDofset(directorDofset &DS) const
{    
	DS.rotation(dofset::tn) = DS.rotation() = DS.rotation(dofset::tn1);
	DS.bodyStepRotationVector().setZero();	
}




void quasiStatic :: incrementDofset(vectorDofset &DS) const
{
	DS.displacement(dofset::tn1) += DS.increment();
	DS.displacement(dofset::tna)  = DS.displacement(dofset::tn1);
	
	extern double global_dt;
	double idt = 1.0/global_dt;
	DS.velocity(dofset::tn1)     += DS.increment() * idt;
	DS.velocity(dofset::tna)      = DS.velocity(dofset::tn1);
}



void quasiStatic :: incrementDofset(scalarDofset &DS) const
{
	DS.displacement(dofset::tn1) += DS.increment();
	DS.displacement(dofset::tna)  = DS.displacement(dofset::tn1);
	
	extern double global_dt;
	double idt = 1.0/global_dt;
	DS.velocity(dofset::tn1)     += DS.increment() * idt;
	DS.velocity(dofset::tna)      = DS.velocity(dofset::tn1);
}



void quasiStatic :: incrementDofset(rotationDofset &DS) const
{
	// rotational dofs
	// rotation qDLambda associated with  rotation vector dTheta: dLambda = exp(dTheta)
	iquaternion qDLambda( DS.increment() );
	
	if (spatialUpdates)
	{
		// left update of the rotation
		DS.rotation(dofset::tn1) = DS.rotation(dofset::tna) = qDLambda * DS.rotation(dofset::tn1);
	}
	else
	{
		// right update of the rotation
		DS.rotation(dofset::tn1) = DS.rotation(dofset::tna) = DS.rotation(dofset::tn1) * qDLambda;
	}
	
	// compute spatial incremental rotation vector from tn to tn+1
	iquaternion qinc = DS.rotation(dofset::tn1)* DS.rotation(dofset::tn).conjugate();
	ivector theta;
	theta.extractFrom(qinc);
	
	// store computed quantities
	DS.spatialStepRotationVector() = theta;
}




// multiplicative update in body rotation 
// d_{n+1}^(k+1) = T_{n+1}^(k+1) E3
//
// T_{n+1}^(k+1) = T_{n+1}^(k) * exp[ dTheta ]
// T_{n+1}^(k+1) = T_{n} * exp[ Theta ]
// 
void quasiStatic :: incrementDofset(directorDofset &DS) const
{
	vector2&	dTheta = DS.increment();
	
#ifndef WITHEIGEN
	iquaternion qdLambda( ivector(dTheta(0), dTheta(1), 0.0) );
#else
	Vector3d vv( dTheta(0), dTheta(1), 0.0 );
	iquaternion qdLambda( AngleAxisd(vv.norm(), vv.normalized() ) );
#endif
	
	// right update of the rotation
	DS.rotation(dofset::tn1) = DS.rotation(dofset::tna) = DS.rotation(dofset::tn1) * qdLambda;
	
	// compute body incremental rotation vector from tn to tn+1
	iquaternion qinc = DS.rotation(dofset::tn1)* DS.rotation(dofset::tn).conjugate();
	ivector Theta;
	Theta.extractFrom(qinc);

	DS.bodyStepRotationVector()[0] = Theta[0];
	DS.bodyStepRotationVector()[1] = Theta[1];
}





void quasiStatic :: info(ostream &of)
{
	integrator :: info(of);
	of << "\n Incremental quasistatic integrator.";
	of << "\n Time is only a parameter to describe the progress of the analysis.";
	of << "\n Rotational updates are done " << (spatialUpdates ? "spatially" : "body");
}







/* this function allocates the space for the linear system of equations
and the rates involved in the solution (none for the static integrator)
*/
void quasiStatic :: initialize(model &model, linearSolver &ls, linearSOE &theSOE, assembler& theAssembler)
{    
	// do general initializations
	integrator :: initialize(model, ls, theSOE, theAssembler);
	theAssembler.setNumberOfTimeDerivatives(0);
	
	
    //  these variables are initialized only once, and not in every step, because
    //  they don't depend on dt
    extern double ftgmass, ftgdamp, ftgstiff;

    ftgstiff = 1.0;
    ftgmass  = 0.0;
    ftgdamp  = 0.0;
	
	DebugMessage("Initializing static integrator");
}




/* this function imposes a dof in a node to have a certain value */
void quasiStatic :: setNodeDof(node &nd, int ndof, double x, double dt)
{	
	if (nd.getNodetype() == _Unode)
	{
		Unode& nnd = dynamic_cast<Unode&>(nd);
		nnd.setDof(ndof, x);
		nnd.getUDS().displacement(dofset::tn1)  = nnd.getUDS().displacement(dofset::tna);
	}
}


// nothing needs to be done, no predictors need to be computed	
void quasiStatic :: startStep(model &m, const double dt)
{}

