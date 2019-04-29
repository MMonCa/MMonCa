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
/* element.cpp
*
* ignacio romero
* may 2000, revised july 2003
* converted to C++ in feb 2006
*
*/

#include "Elements/element.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <limits>
#include <stdexcept>

#include <cassert>
#include <cctype>
#include <cstdlib>
#include <cstdarg>
#include <cmath>
#include <vector>
#include <limits>

#include "Elements/eltype.h"
#include "Elements/Interpolation/interpolation.h"
#include "Elements/evalspot.h"

#include "Model/Node/node.h"
#include "Materials/material.h"

#include "Analysis/assembler.h"
#include "Analysis/analysis.h"
#include "Analysis/Integrators/integrator.h"
#include "Analysis/Integrators/quasistatic.h"

#include "Geometry/quadpoint.h"
#include "General/idata.h"

#include "Io/logger.h"
#include "Io/message.h"

#include "Math/feliksmath.h"


using  namespace std;
using namespace blue;


element :: element(const int labelc, const int type, const feliks::topology::cell* c, const int nn, node **ndlist) :
	active(true),
	eigmax(0.0),
	time(0.0),
	etype(0),
	typelabel(type),
	label(labelc),
	lumpedMass(false),
	hasConstrDofs(false),
	nfaces(0),
	mass(0.0),
	theMaterialDescription(NULL),
	parentSubmesh(0),
	theEnergies(),
    theCell(c)
{
    nodes.reserve(nn);
    for (int a=0; a<nn; a++) nodes.push_back(ndlist[a]);
}



element :: element(const int l, const eltype &t, const feliks::topology::cell* c, const int nn, node **ndlist)
:	label(l),
	eigmax(0.0),
	time(0.0),
	etype(&t),
	active(true),
	lumpedMass(false),
	hasConstrDofs(false),
	nfaces(0),
	mass(0.0),
	theMaterialDescription(NULL),
    parentSubmesh(0),
    theCell(c)
{
    nodes.reserve(nn);
    for (int a=0; a<nn; a++) nodes.push_back(ndlist[a]);
	typelabel = t.label();	
}


element :: element(const int labelx) 
:	label(labelx),
	eigmax(0),
	time(0.0),
	etype(0),
	nodes(0),
	active(true),
	lumpedMass(false),
	hasConstrDofs(false),
	nfaces(0),
	mass(0.0),
	theMaterialDescription(NULL),
	parentSubmesh(0),
	typelabel(0),
    theCell(0)
{
    
}




/* deallocates the memory asigned for the element. The nodes stored in the list
* c->nodes, must be deallocated elsewhere. This is because we might want to preserve
* the nodes for other elements.
*/
element :: ~element()
{	
    neighbors.clear();
	nodes.clear();
}




void element :: addConnectivity(element &e2)
{
    this->neighbors.push_back(&e2);
    e2.neighbors.push_back(this);
}




// this function finds the eigenvalues of the generalized problem K v = w^2 M v, where M is the
// lumped mass. To do so, it inverts the sqrt mass matrix and multiplies K by it, obtaining
// a classical EV problem
void element ::	 approximatedEigenvalues() const
{
	matrix K;
	K.resize( getNDofs(), getNDofs() );
	K.setZero();
	element* noconst = const_cast<element*>(this);
	staticTangent(*noconst, K);
	
	// scale. Premultiply and postmultiply tb by M^half
	size_t ndofs = nodes[0]->ndofs;
	for (int a=0; a<getNNodes(); a++) 
		for (int b=0; b<getNNodes(); b++) 
			for (int i=0; i<ndofs; i++)
				for (int j=0; j<ndofs; j++)
					K.data[a*ndofs+i][b*ndofs+j] *= sqrt(nodes[a]->getInvMass()) * sqrt(nodes[b]->getInvMass());


	
    matrix evectors(getNDofs(), getNDofs());
    longvector evaluesR(getNDofs());
    longvector evaluesI(getNDofs());  // imaginary part.
	
	K.eigendata(evectors, evaluesR, evaluesI);
	
	vector<double> eval;
	for (int a=0; a<getNDofs(); a++) eval.push_back(evaluesR.data[a]);
	sort(eval.begin(), eval.end());
    
    Message("\n Generalized eigenvalues for element %d and corresponding time step sizes", getLabel());
	Message("\n Nodal masses are not recomputed, so they come from the assembly of all neighbor elements");
	Message("\n     eigenvalues w^2            w            dt=2/w:");
    for (int a= (int)getNDofs()-1; a>=0; a--)
    {
        Message("\n     %2d : %+10.4e      %+10.4e      %+10.4e", a+1, 
			eval[a], sqrt(eval[a]), 2.0/sqrt(eval[a]) );
    }
	
	// Recompute the nodal masses, so that only the element mass is added
	vector<double> oldMasses;
	for (size_t k=0; k<getNNodes(); k++) oldMasses.push_back(nodes[k]->getMass());
	
	for (size_t k=0; k<getNNodes(); k++) nodes[k]->resetMass();
	noconst->lumpMassToNodes();

	// compute tangent again
	K.setZero();
	staticTangent(*noconst, K);
	
	// scale. Premultiply and postmultiply tb by M^half
	for (int a=0; a<getNNodes(); a++) 
		for (int b=0; b<getNNodes(); b++) 
			for (int i=0; i<ndofs; i++)
				for (int j=0; j<ndofs; j++)
					K.data[a*ndofs+i][b*ndofs+j] *= sqrt(nodes[a]->getInvMass()) * sqrt( nodes[b]->getInvMass() );


	K.eigendata(evectors, evaluesR, evaluesI);
	eval.clear();
	for (int a=0; a<getNDofs(); a++) eval.push_back(evaluesR.data[a]);
	sort(eval.begin(), eval.end());
    
    Message("\n\n Generalized eigenvalues for element %d and corresponding time step sizes", getLabel());
	Message("\n Nodal masses are recomputed, and only the element contribution is added");
	Message("\n     eigenvalues w^2            w            dt=2/w:");
    for (int a=(int)getNDofs()-1; a>=0; a--)
    {
        Message("\n     %2d : %+10.4e      %+10.4e      %+10.4e", a+1, 
			eval[a], sqrt(eval[a]), 2.0/sqrt(eval[a]) );
    }

	
	double dta(noconst->maxEigenvalue());
	Message("\n\n Estimated max eigenvalue upper bound w = %10.4e", dta);
	Message("\n Time step size estimate used by feliks dt = %10.4e", 2.0/dta);

	// reconstruct masses
	for (size_t a=0; a< getNNodes(); a++) 
	{
		nodes[a]->resetMass();
		nodes[a]->incrementMass( oldMasses[a] );
	}
	 
}






/* this function calls the check function of each formulation. If none is defined,
 by default it returns true, which means there is no error
*/
bool element :: check()
{
	return true;
}



// shortest distance between nodes
double element :: characteristicDim()
{
    double h=std::numeric_limits<double>::max();
    
    vector<node*>::iterator  it1 = nodes.begin();
    while ( it1 != nodes.end() )
    {
        vector<node*>::iterator  it2 = it1+1;
        while (it2 != nodes.end() )
        {
            double d = ( (*it1)->coordinates() - (*it2)->coordinates() ).norm();
            h = std::min<double>( h , d);
            ++it2;
        }
        ++it1;
    }
    
    return h;
}






/*  this function builds the '*faces' array. 
*	The "face" concept here is tricky. We want to compute faces for contact and thus a "face" is
*	not the boundary of the element, but depends on the surrounding space. The same geometry 
*	(a quad, for example) can have different faces depending on the space it lives on. In 2d, 
*	its faces are its 4 edges. In 3d, it has two faces, the quad itself but from with the two
*	opposite orientations.  It has to be done like this for the contact. This function must
*	be programmed together with ElementNodeLabelsInFace(), which returns the local node numbers
*	for each face.
* 
*  In each entry of the *faces array there is an idata which
*  contains the labels of the nodes in that face. The numbering order of the faces is important
*  because we will later need to recognize the faces of an element.
*  We use a double pointer for faces because this pointer is defined before calling
*  the function, and we must only modify its content, not its value.
*   
*/
int element :: defineFaces(size_t *nfaces, idata **faces)
{
    geometryT geo;
    size_t    npf=0, k; //npf = nodes per face
	const	int ndm(3);

	geo   = etype->geometry();
    
    switch (1000*ndm + 100*geo + getNNodes())
    {
        // linear elements with 2 nodes in 1 dimension
        case 1102:
            *nfaces = 2;
            npf     = 1;
            break;
			
		// linear elements in 2d
		case 2102:
			*nfaces = 2;
			npf     = 2;
			break;
			
		// 2d triangle and quad
        case 2203:
		case 2204:
            *nfaces = getNNodes();
            npf     = 2;
            break;
			
		//3d triangle and 3d quad
		case 3203:
		case 3204:
			*nfaces = 1;
			npf     = getNNodes();
			break;
			
		// 3d tet
        case 3304:
            *nfaces = 4;
            npf     = 3;
            break;
			
		// 3d hex
        case 3308:
            *nfaces = 6;
            npf     = 4;
            break;
            
        default:
            WarningMessage("ElementDefineFaces not programmed for element with %d nodes in %dd", getNNodes(), ndm);
            *nfaces = 0;
            break;
    }
    
    if (*nfaces >0)
    {
        // allocate space for all the idatas that must hold the node labels
        *faces = (idataT **) malloc( (*nfaces) * sizeof(idata));
        
        // loop through all the faces, filling the data of node labels
        for (k=0; k< *nfaces; k++)
        {
            (*faces)[k] = NewIdata((int) npf);
            getNodeLabelsInFace(k, (*faces)[k]);
        }
    }
    
    return 1;
}




/* computation of the numeric tangent from the centered difference quotient
 */
bool  element :: computeNumericDDEnergy(matrix &tg)
{
    double   inc = 1e-6;
	singleElementAssembler eas(*this);
    double   alpha = 1.0;
    
    // loop through all the dofs computing the quotient
	eas.getTasks().reset();
	eas.emat.setZero();
	eas.getTasks().static_residual = true;
    
    for (int b=0; b<getNNodes(); b++)
    {
		if (getNode(b).getNodetype() == _Unode  || 
            getNode(b).getNodetype() == _UPnode ||
            getNode(b).getNodetype() == _UDnode )
		{
			Unode&   nd(  dynamic_cast<Unode&>(getNode(b)) );
			ivector& ub_eq = nd.getUDS().displacement();
			ivector& ub_n1 = nd.getUDS().displacement(dofset::tn1);
			
			int i=0;
			while (i<3 && !nd.getUDS().isConstrained(i))
			{
				double original_eq = ub_eq(i), original_n1 = ub_n1(i);
				
				ub_eq(i) += alpha*inc;
				ub_n1(i) +=       inc;
				eas.eres.setZero();
				contributeToGenericIntegral(eas);
				longvector rplus = eas.eres;
				
				ub_eq(i) -= 2.0*alpha*inc;
				ub_n1(i) -= 2.0*inc;
				eas.eres.setZero();
				contributeToGenericIntegral(eas);
				longvector rminus = eas.eres;
                
				ub_eq(i) = original_eq;
				ub_n1(i) = original_n1; 
				const int column = b*getNDofsPerNode()+i;
				for (int k=0; k< getNDofs(); k++)
				{
					const double Tij = -(rplus(k)-rminus(k))/(2.0*inc);
					eas.emat(k, column) = Tij;
				}
				++i;
			}
		}
        
		if ( getNode(b).getNodetype() == _UPnode )
		{
			UPnode&   nd(  dynamic_cast<UPnode&>(getNode(b)) );
			double*   ub_eq = &nd.getPDS().displacement();
			double*	  ub_n1 = &nd.getPDS().displacement(dofset::tn1);
			
			
			if ( !nd.getPDS().isConstrained())
			{
				double original_eq = *ub_eq;
				double original_n1 = *ub_n1;
				
				*ub_eq += alpha*inc;
				*ub_n1 +=       inc;
				eas.eres.setZero();
				contributeToGenericIntegral(eas);
				longvector rplus = eas.eres;
				
				
				*ub_eq -= 2.0*alpha*inc;
				*ub_n1 -= 2.0*inc;
				eas.eres.setZero();
				contributeToGenericIntegral(eas);
				longvector rminus = eas.eres;
				
				
				*ub_eq = original_eq;
				*ub_n1 = original_n1; 
				const int column = b*getNDofsPerNode()+3;
				for (int k=0; k< getNDofs(); k++)
				{
					const double Tij = -(rplus(k)-rminus(k))/(2.0*inc);
					std::cout<< "T(" << k << ", " << column << ") = " << Tij <<std::endl;
					eas.emat(k, column) = Tij;
				}
			}	
			
		}
        
        
		if ( getNode(b).getNodetype() == _UDnode )
		{
            quasiStatic theIntegrator;
            inc = 1e-3;
            
			UDnode&   nd(  dynamic_cast<UDnode&>(getNode(b)) );
			iquaternion r_eq = nd.getDDS().rotation();
			iquaternion r_n1 = nd.getDDS().rotation(dofset::tn1);
            
			int i=0;
			while (i<2 && !nd.getDDS().isConstrained(i))
			{
				nd.getDDS().increment().setZero();
				nd.getDDS().increment()[i] = inc;
				theIntegrator.incrementDofset( nd.getDDS() );
				eas.eres.setZero();
				contributeToGenericIntegral(eas);
				longvector rplus = eas.eres;
				nd.getDDS().rotation(dofset::tn1) = r_n1;
				nd.getDDS().rotation(dofset::tna) = r_eq;
                
				nd.getDDS().increment().setZero();
				nd.getDDS().increment()[i] = -inc;
				theIntegrator.incrementDofset( nd.getDDS() );
				eas.eres.setZero();
				contributeToGenericIntegral(eas);
				longvector rminus = eas.eres;
				nd.getDDS().rotation(dofset::tn1) = r_n1;
				nd.getDDS().rotation(dofset::tna) = r_eq;
				nd.getDDS().increment().setZero();
				
				const int column = b*getNDofsPerNode()+3+i;
				for (int k=0; k< getNDofs(); k++)
				{
					const double Tij = -(rplus(k)-rminus(k))/(2.0*inc);
					eas.emat(k, column) = Tij;
				}
				++i;
			}
		}
	}
    
	tg = eas.emat;
    
    /*
     // store the status of the node before the numerical derivation. In the
     // future, this should also include the rates
     dofs = nd.getDofs();
     
     Dcopy(dofs, Znddofs->data, ndofs);
     if (nd.rotationType() > ROT_NONE) Zrotation = nd.rotationMatrix();
     
     for (int i=0; i<ndofs; i++)
     {
     double original = nd.
     
     // this is the standard way to compute the numeric tangent. However, for nonlinear
     // updates (rods and shells, we need something more sophisticated
     double Zu = nd.getDof(i);
     
     // compute forward residual
     nd.setDof(i,  Zu + inc);
     nd.updateCurrentCoordinates();
     ZeroVector(rplus);
     residual(rplus, 'a');
     
     // compute backward residual
     nd.setDof(i,  Zu - inc);
     ZeroVector(rminus);
     nd.updateCurrentCoordinates();
     residual(rminus, 'a');
     
     // reset the original (unperturbed) dof
     nd.setDof(i, Zu);
     nd.updateCurrentCoordinates();
     
     // the following version should work even for non-additive nodal updates because the update is not done here.
     // Also, it can/will include the updates in the rates
     
     // forward residual
     ZeroVector(rplus);
     nd.setDeltaDofsToZero();
     nd.setDofIncrement(i, inc);
     theIntegrator.incrementNodeSolution(nd);
     residual(rplus, 'a');
     nd.setDofs(Znddofs);
     if (nd.rotationType() > ROT_NONE) nd.setRotation(Zrotation);
     
     // backward residual
     ZeroVector(rminus);			 
     nd.setDeltaDofsToZero();
     nd.setDofIncrement(i, -inc);
     theIntegrator.incrementNodeSolution(nd);
     residual(rminus, 'a');
     nd.setDofs(Znddofs);
     if (nd.rotationType() > ROT_NONE) nd.setRotation(Zrotation);
     
     for (int k=0; k<getNNodes()*ndofs; k++)
     tg->data[k][a*ndofs+i] = -(rplus->data[k]-rminus->data[k])/(2.0*inc);
     }	
     }
     
     FreeVector(Znddofs);
     FreeMatrix(Zedofs);
     FreeVector(rplus);
     FreeVector(rminus);
     */
	return 1;
}





bool element ::	contributeToGenericIntegral(assembler& theAssembler)
{
	return sharedFunction(theAssembler);	
}



void element :: contributeToResidual(assembler& theAssembler)
{
	// assemble standard contribution
	sharedFunction(theAssembler);
	
	// add modifications due to imposed bc
	extern int global_iteration;
	if (global_iteration == 0 && hasConstrainedDofs())
	{
		BCassembler bca(theAssembler);
		sharedFunction(bca);
	}
}




void element :: contributeToResidualTangent(assembler& theAssembler)
{	
	sharedFunction(theAssembler);
	
	extern int global_iteration;
	if (global_iteration == 0 && hasConstrainedDofs())
	{
		matrix dummytangent(0,0);
		BCassembler bca(theAssembler);
		sharedFunction(bca);		
	}	
}




// this function is required for the sort in the following functions
bool element :: epLess(const element *ep1, const element *ep2) 
{
	return *ep1 < *ep2;
}





/* assemble the element residual fext-fint and store it in the nodal
variables. The fint vector only contains the part that must be evaluated
in an explicit method, that is, it does not contain the terms of the
higher time derivate. In case of solid dynamics, it does not contain
inertial terms

// temporary defintion of the function. Since function is virtual, it should provide
// the implementation for all these elements that don't have the function yet.
// the way it does it is by using the 'residual' function. This is completely
// equivalent to what we want to do, but much slower. We should replace this implementation
// with more efficient ones in each formulation

bool element :: explicitForcesToNodes()
{
    //  used degrees of freedom per node of the element. Careful: nodes might
	//have more degrees of freedom, because they belong to other elements. At
	//least they must have this number of dofs: 
	int      ndofn( etype->getDofsPerNode() );
	
    // get all the element forces = fext - fint
	matrix dummy=NULL;
	elmtTasks tasks;
	tasks.static_residual = true;
	tasks.dynamic_residual = false;
	longvector res=NULL;
	res = NewVector( getNDofs() );
	sharedFunction(res, dummy, tasks);
	
	longvector nodef=NULL;
    for (int i=0; i< getNNodes(); i++)
    {
        // extract the forces corresponding to node i 
        nodef = SubVector(res, i*ndofn, (i+1)*ndofn-1);
		
        // and insert them into the node structure
		node &nd = getNode(i);
        nd.assembleForce(nodef);
		
        FreeVector(nodef);
    }
	
    FreeVector(res);
	return true;
}

*/



bool element :: explicitForcesToNodes()
{
	logger::mainlog << endl << "Element " << getLabel() 
        << " is not prepared for explicit analysis (type = " << etype->name() << ")" << flush;
	return false;
}





/* this retrieves a list that includes all the dofs number in the element,
* without any attention to which node they belong.
*/
void  element :: getDofNumbers(vector<int>& theDofs)
{
    //  these are the dofs/node that the element actually requires
    int nddofs = etype->getDofsPerNode();
    	
    // now we fill up the idata with the degrees of freedom
    theDofs.clear();
    theDofs.reserve(nddofs*getNNodes());
    for (int k=0; k<getNNodes(); k++)
    {
        for (int j=0; j<nddofs; j++)
			theDofs.push_back(nodes[k]->getID(j));
    }
}




std::set<evalspot*> element :: getEvalspots() const
{
    std::set<evalspot*> junk;
    return junk;
}





/*  This function fills up the nodes vector with the pointers to the nodes of the
face 'f' of element 'e'. 
*/
void element :: getNodesInFace(const int face, std::vector<node*> &faceNodes)
{
	for (int k=0; k< getNNodesPerFace(); k++)
		faceNodes.push_back(&getNodeInFace(face, k));
}




node& element :: getNodeInFace(const int face, const int local) const
{
	const int ndm(3);
	node *nd(0);
	
    switch (100000*ndm + 10000*etype->geometry() + 100*getNNodes() + 10*face + local)
    {
        // 1d/2d/3d linear elements. Faces are points
		case 110200:
		case 210200:
		case 310200: nd = nodes[0]; break;
			
		case 110210:
		case 210210:
		case 310210: nd = nodes[1]; break;
			
			
		// 2d triangle. Faces are segments
        case 220300: nd = nodes[0]; break;
		case 220301: nd = nodes[1]; break;

		case 220310: nd = nodes[1]; break;
		case 220311: nd = nodes[2]; break;

		case 220320: nd = nodes[2]; break;
		case 220321: nd = nodes[0]; break;			
			
		// 3d triangle. Faces are triangle
		case 320300: nd = nodes[0]; break;
		case 320301: nd = nodes[1]; break;
		case 320302: nd = nodes[2]; break;
			
            
		// 2d quad. Faces are segments
        case 220400: nd = nodes[0]; break;
        case 220401: nd = nodes[1]; break;
			
        case 220410: nd = nodes[1]; break;
        case 220411: nd = nodes[2]; break;

		case 220420: nd = nodes[2]; break;
        case 220421: nd = nodes[3]; break;

		case 220430: nd = nodes[3]; break;
        case 220431: nd = nodes[0]; break;			
			
			
		// 3d quad. Faces are quads
		case 320400: nd = nodes[0]; break;
		case 320401: nd = nodes[1]; break;
		case 320402: nd = nodes[2]; break;
		case 320403: nd = nodes[3]; break;
			
			
		// tet. Faces are triangles
		case 330400: nd = nodes[0]; break;
		case 330401: nd = nodes[1]; break;
		case 330402: nd = nodes[3]; break;
			
		case 330410: nd = nodes[1]; break;
		case 330411: nd = nodes[2]; break;
		case 330412: nd = nodes[3]; break;

		case 330420: nd = nodes[2]; break;
		case 330421: nd = nodes[0]; break;
		case 330422: nd = nodes[3]; break;

		case 330430: nd = nodes[0]; break;
		case 330431: nd = nodes[2]; break;
		case 330432: nd = nodes[1]; break;
			
			
		// hex. Faces are quads
		case 330800: nd = nodes[0]; break;
		case 330801: nd = nodes[3]; break;
		case 330802: nd = nodes[2]; break;
		case 330803: nd = nodes[1]; break;

		case 330810: nd = nodes[0]; break;
		case 330811: nd = nodes[1]; break;
		case 330812: nd = nodes[5]; break;
		case 330813: nd = nodes[4]; break;

		case 330820: nd = nodes[1]; break;
		case 330821: nd = nodes[2]; break;
		case 330822: nd = nodes[6]; break;
		case 330823: nd = nodes[5]; break;

		case 330830: nd = nodes[2]; break;
		case 330831: nd = nodes[3]; break;
		case 330832: nd = nodes[7]; break;
		case 330833: nd = nodes[6]; break;

		case 330840: nd = nodes[0]; break;
		case 330841: nd = nodes[4]; break;
		case 330842: nd = nodes[7]; break;
		case 330843: nd = nodes[3]; break;

		case 330850: nd = nodes[4]; break;
		case 330851: nd = nodes[5]; break;
		case 330852: nd = nodes[6]; break;
		case 330853: nd = nodes[7]; break;
						
        default:
            WarningMessage("ElementNodesInFace not programmed for element with %d nodes in %dd", getNNodes(), ndm);
            break;
    }
	return *nd;
}





void  element :: getNodeLabels(vector<int>& labels) const
{
    labels.clear();
    for (int k=0; k< getNNodes(); k++) 
        labels.push_back(nodes[k]->getLabel());	
}





// put all the labels of the nodes on an element face inside the idata
// faces are oriented according to right hand rule, with outward normal
// the labels must be a non-null idata
void element :: getNodeLabelsInFace(int face, idata &labels) const
{
    geometryT geo;
	const int ndm(3);
	
	IdataEmpty(labels);
    
    geo = etype->geometry();
    
    switch (10000*ndm + 1000*geo + 10*getNNodes() + face)
    {
        // 1d/2d/3d linear elements. Faces are points
        case 1020: AppendToIdata(labels, nodes[0]->getLabel() ); break;
        case 1021: AppendToIdata(labels, nodes[1]->getLabel() ); break;
            
		// 2d triangle. Faces are segments
        case 22030: IdataAppend(labels, 2, nodes[0]->getLabel(), nodes[1]->getLabel() ); break;
        case 22031: IdataAppend(labels, 2, nodes[1]->getLabel(), nodes[2]->getLabel() ); break;
        case 22032: IdataAppend(labels, 2, nodes[2]->getLabel(), nodes[0]->getLabel() ); break;
			
		// 2d quad. Faces are segments
        case 22040: IdataAppend(labels, 2, nodes[0]->getLabel(), nodes[1]->getLabel() ); break;
        case 22041: IdataAppend(labels, 2, nodes[1]->getLabel(), nodes[2]->getLabel() ); break;
        case 22042: IdataAppend(labels, 2, nodes[2]->getLabel(), nodes[3]->getLabel() ); break;
        case 22043: IdataAppend(labels, 2, nodes[3]->getLabel(), nodes[0]->getLabel() ); break;
			
		// 3d triangle. Faces are the triangle itself
		case 32030: IdataAppend(labels, 3, nodes[0]->getLabel(), nodes[1]->getLabel(), nodes[2]->getLabel() ); break;
			
		// 3d quad. Faces are the quad itself
        case 32040: IdataAppend(labels, 4, nodes[0]->getLabel(), nodes[1]->getLabel(), nodes[2]->getLabel(), nodes[3]->getLabel() ); break;
			
		// tet. Faces are triangles
        case 33040: IdataAppend(labels, 3, nodes[0]->getLabel(), nodes[1]->getLabel(), nodes[3]->getLabel() ); break;
        case 33041: IdataAppend(labels, 3, nodes[1]->getLabel(), nodes[2]->getLabel(), nodes[3]->getLabel() ); break;
        case 33042: IdataAppend(labels, 3, nodes[2]->getLabel(), nodes[0]->getLabel(), nodes[3]->getLabel() ); break;
        case 33043: IdataAppend(labels, 3, nodes[0]->getLabel(), nodes[2]->getLabel(), nodes[1]->getLabel() ); break;
			
		// hex. Faces are quads
        case 33080: IdataAppend(labels, 4, nodes[0]->getLabel(), nodes[3]->getLabel(), nodes[2]->getLabel(), nodes[1]->getLabel() ); break;
        case 33081: IdataAppend(labels, 4, nodes[0]->getLabel(), nodes[1]->getLabel(), nodes[5]->getLabel(), nodes[4]->getLabel() ); break;
        case 33082: IdataAppend(labels, 4, nodes[1]->getLabel(), nodes[2]->getLabel(), nodes[6]->getLabel(), nodes[5]->getLabel() ); break;
        case 33083: IdataAppend(labels, 4, nodes[2]->getLabel(), nodes[3]->getLabel(), nodes[7]->getLabel(), nodes[6]->getLabel() ); break;
        case 33084: IdataAppend(labels, 4, nodes[3]->getLabel(), nodes[0]->getLabel(), nodes[4]->getLabel(), nodes[7]->getLabel() ); break;
        case 33085: IdataAppend(labels, 4, nodes[4]->getLabel(), nodes[5]->getLabel(), nodes[6]->getLabel(), nodes[7]->getLabel() ); break;
			
        default:
            WarningMessage("ElementNodeLabelsInFace not programmed for element with %d nodes in %dd", getNNodes(), ndm);
            break;
    }
}






/* this is a tricky question because one element might have nodes with different
* number of degrees of freedom: for example, a node that connects to a beam (6 dof)
* and a node that connects to a solid (3 dof). What this function answers is
* the number of dofs associated with the formulation implemented in the element.
* This is as if all the nodes had the same number of degrees of
* freedom. This number is given in the eltype definition.
*/
size_t element :: getNDofs() const
{
    return (getNNodes() * etype->getDofsPerNode() );
}




size_t element ::	getNDofsPerNode() const
{
	return etype->getDofsPerNode();
}



/* compute the number of faces (hypersurfaces) in the boundary of an element
*/
size_t element :: nFaces() const
{
    size_t nf=0;
	const int ndm(3);
	
    switch (1000*ndm + 100*etype->geometry() + getNNodes())
    {
        // linear elements with 2 nodes in 1 dimension
        case 1102: 
		case 3102: nf = 2; break;
			
		// linear elements in 2d
		case 2102: nf = 2; 	break;
			
		// 2d triangle and quad
        case 2203:
		case 2204: nf = getNNodes(); break;
			
		//3d triangle and 3d quad
		case 3203:
		case 3204: nf = 1; break;
			
		// 3d tet
        case 3304: nf = 4; break;
			
		// 3d hex 
        case 3308: nf = 6; break;
            
        default:
            WarningMessage("getNFaces not programmed for element with %d nodes in %dd", getNNodes(), ndm);
    }	
    return nf;
}




element* element :: getNeighbor(int n) const
{
    assert( n < nNeighbors());
    return neighbors[n];
}



int element :: getNNodesPerFace() const
{
	geometryT geo;
    int       npf=0;
	int ndm(3);

    geo = etype->geometry();
    
    switch (1000*ndm + 100*geo + getNNodes())
    {
        // linear elements with 2 nodes
		case 1102:  case 2102 : case 3102: npf = 1; break;
            
		// 2d/3d triangle and quad
        case 2203: case 2204: npf = 2; break;
			
		//3d triangle, 3d quad
		case 3203: 
		case 3204: npf = getNNodes(); break;

		// 3d tet
        case 3304: npf = 3; break;
			
		// 3d hex
        case 3308: npf = 4; break;
            
        default:
            WarningMessage("ElementNNodesPerFace not programmed for element with %d nodes in %dd", getNNodes(), ndm);
            break;
    }
	return npf;
}




/* we obtain the number of quadrature points looking at the number of nodes and the
geometry type of the element. That is because the spatial dimension of the problem
is now what matters. If an element has the geometry of a curve (1), the number
of quadrature points is independent of the fact that the element is in a 1,2, or 3
dimensional space */
int element :: nQuadraturePoints() const
{
    return numberOfQuadraturePoints(getNNodes(), etype->geometry() );
}






/* this function "guesses" the maximum degree of complete polynomials in the element
interpolation. This is just a guess because then the particular element may be
subintegrated. This guess can be obtained solely from the dimension and # of nodes.
This routine is used for the Superconvergence Patch Recovery error estimation or
just for the stress smoothing.
*/
int element :: getPolynomialDegree()
{
    int d;
	int ndm(3);

	
    switch ( 100*etype->geometry() + getNNodes() )
    {
        case 102: d = 1; break;	/* linear 1d element */
        case 103: d = 2; break;	/* quadratic 1d element */
        case 203: d = 1; break;	/* cst  */
        case 206: d = 2; break; /* 6 noded quadratic triangle */
        case 204: d = 1; break;	/* quad */
        case 208: d = 2; break;	/* serendipity */
        case 209: d = 2; break;	/* quadratic quad */
        case 304: d = 1; break;	/* tet */
        case 310: d = 1; break;	/* quadratic tet */
        case 308: d = 1; break;	/* brick */
        default:
            Message("Error in ElementPolynomialDegree. Element with %d nodes in ndm=%d not defined",
                    getNNodes(), ndm);
            d = -1;
    }
    return d;
}





/* fills the element data with information from the geometry and input file
already stored.
*/
bool element :: initialize()
{
	//  if the element has geometric entity (not points, nor contact elements)
	//  compute number of faces and allocate space for as many neighbors as possible.
	if (etype->geometry() > POINT)
	{
		nfaces   = nFaces();
        neighbors.reserve(nfaces);
	}
	
	
	// links to material type and creates material points
	int matlabel  = etype->materialLabel();
	
	if (matlabel > 0) 	
	{
		theMaterialDescription = &(material::getMaterial(matlabel));
	}
	 
	
    // link nodes to element and update number of dofs in each node
    int       dofsnode  = etype->getDofsPerNode();
    //rotationT rt        = etype->getRotationType();
	
    for (int k=0; k<getNNodes(); k++)
    {
        nodes[k]->setNdofs(dofsnode);
		
		/*
		// if there are any rotations,  allocate space and set to defaults
		ivector E3(0.0, 0.0, 1.0);
		iquaternion q;
        if      (rt == ROT_SHELL)		nodes[k]->initializeRotationalDofs(E3);
		else if (rt == ROT_ROD_EXPMAP)	nodes[k]->initializeRotationalDofs(q);
		 */
    }
	
	return true;
}



bool element :: integrateDampingMatrix(assembler& as) 
{
    cout << "integrateDampingMatrix in element not programmed yet";
    return false;
}




bool element :: integrateMassMatrix(assembler& as)
{
	return sharedFunction(as);	
}



// this function interpolates the acceleration to a point as dictated by the shp
void element :: interpolateAcceleration(const FEshapefunbuilder& shp, ivector &acc, const dofset::evaluation_time when)
{
	acc.setZero();
	
	for (int a=0; a< getNNodes(); a++) 
		acc += shp[a].value * dynamic_cast<Unode&>( getNode(a) ).getUDS().acceleration(when);
}


/*

// this function interpolates the spatial angular velocity in problems that have rotational
// degress of freedom (this is, 6 dofs)
// It does not work for plane problems !!
void element :: interpolateAngularVelocity(const FEshapefunbuilder& shp, ivector& omega)
{
	omega.setZero();
	for (int a=0; a<getNNodes(); a++) omega += shp[a].value * nodes[a]->angularVelocity();
}

*/

void element :: interpolateBodyAlpha(const FEshapefunbuilder& shp, ivector& Alpha)
{
}


// this function interpolates the body angular velocity in problems that have rotational
// degress of freedom (this is, 6 dofs)
// It does not work for plane problems !!
void element :: interpolateBodyOmega(const FEshapefunbuilder& shp, ivector& Omega)
{
}



// this function interpolates the body rotation vector that takes from Lambda_n
// to Lambda_n+1
void element :: interpolateBodyRotationVector(const FEshapefunbuilder& shp, ivector& Theta)
{
	Theta.setZero();
	//ROT for (int a=0; a< getNNodes(); a++) Theta += shp[a].value * nodes[a]->bodyRotationVector();
}




double element :: interpolatePressure(const FEshapefunbuilder& shp, const dofset::evaluation_time when) const
{
	if ( getNode(0).getNodetype() != _UPnode)
		throw std::runtime_error("Trying to interpolate pressure in a nonpressure element");
	
	double p(0.0);

	if (when != dofset::t0)
		for (int a=0; a<getNNodes(); a++) 
			p += dynamic_cast<UPnode&>(getNode(a)).getPDS().displacement(when) * shp[a].value;


	return p;
}



double element :: interpolatePressureRate(const FEshapefunbuilder& shp, const dofset::evaluation_time when)
{
	if ( getNode(0).getNodetype() != _UPnode)
		throw std::runtime_error("Trying to interpolate pressure in a nonpressure element");
	
	double p(0.0);

	for (int a=0; a<getNNodes(); a++) 
			p += dynamic_cast<UPnode&>(getNode(a)).getPDS().velocity(when) * shp[a].value;
	
	return p;
}


double element :: interpolatePressureRate2(const FEshapefunbuilder& shp, const dofset::evaluation_time when)
{
	if ( getNode(0).getNodetype() != _UPnode)
		throw std::runtime_error("Trying to interpolate pressure in a nonpressure element");
	
	double p(0.0);
	
	for (int a=0; a<getNNodes(); a++) 
		p += dynamic_cast<UPnode&>(getNode(a)).getPDS().acceleration(when) * shp[a].value;
	
	return p;
}



void element :: interpolatePosition	(const FEshapefunbuilder& shp, 
                                     ivector &pos, 
                                     const dofset::evaluation_time when) const
{
	pos.setZero();
	for (int a=0; a<getNNodes(); a++) pos += shp[a].value * nodes[a]->coordinates(when);
}






void element :: interpolateReferencePosition(const FEshapefunbuilder& shp, ivector &pos)
{
	pos.setZero();
	for (int a=0; a<getNNodes(); a++) pos += shp[a].value * nodes[a]->getReferenceCoordinates();
}




// interpolates the vector that takes from Lambda_{n+1}^k to Lambda_{n+1}^k+1, i.e., the last iterate
//
// June 24 2009, I have changed the way the incremental vectors are recuperated. In the old feliks,
// when the iteration is 0, there are some dU, but one only gets zero. Now you get what it is stored
//
// May 2011: to make displacement and rotation updates similar, now in the first iteration the interpolated
// rotation incremetn is always zero. The function trickyRotationIncrementVector is the key
void element :: interpolateRotationIncrementVector(const FEshapefunbuilder& shp, ivector& dtheta, ivector& dthetapr)
{
}




void element :: interpolateRotationVector(const FEshapefunbuilder& shp, ivector& theta)
{
	theta.setZero();
    std::cout << "\n Function interpolateRotationVector not implemented";
	//for (int a=0; a< getNNodes(); a++) theta += shp[a].value * nodes[a]->rotationVector();
}






// this function interpolates the velocity to a point as dictated by the shp
void element :: interpolateVelocity(const FEshapefunbuilder& shp, 
                                    ivector &vel, 
                                    const dofset::evaluation_time when) const
{	
	vel.setZero();	
	for (int a=0; a<getNNodes(); a++) vel += shp[a].value * nodes[a]->velocity(when);
}





// we provide a default implementation of lumpDampingToNodes so that elements with no
// mass don't have to implement it
bool element :: lumpDampingToNodes()
{
	return true;
}


// we provide a default implementation of lumpMassToNodes so that elements with no
// mass don't have to implement it
bool element :: lumpMassToNodes()
{
	return true;
}


// figure out if any of the nodes of the element has constrained
// dofs, and set the corresponding flag
// This is used to speed up the residual assembly.
void element :: markConstrained()
{
	hasConstrDofs = false;
	for (int k=0; k<getNNodes(); k++)
	{
		if (nodes[k]->getNConstrainedDofs() > 0)
		{
			hasConstrDofs = true;
			break;
		}
	}
}




void  element :: print()
{
	string name;
	name = etype->name();
	
	Message("\nElement   %4d, ", label);
    Message(  "Type label : %d, ",  typelabel);
    Message(  "Type       : %s, ",  name.c_str()   );
	if (nfaces > 0) Message("\n\t     Number of faces : %4d", nfaces);
	Message("\n\t     Neighbors(%d)", nNeighbors() );
	for (int k=0; k<nNeighbors(); k++) Message(" %3d,", neighbors[k]->getLabel() );
	Message(  "  Nodes : ");
	
	for (int k=0; k<getNNodes(); k++)
	{
		node *nd  = nodes[k];
		if  (nd  == NULL)  Message("UNDEF\t");
		else               Message(" %4d", nd->getLabel() );
	}
}







/* dumps the gradients of an element to file, with no labels.
This function is for the postprocessing, or the loggers
*/
void  element :: printGradients(ostream& of) const
{ 
	cout << "Element::printGradients is not implemented" << endl;
}



void element :: printTgEigenv()
{
    int n, a;
    double im;
	
    // total number of dofs in the element
    n = getNNodes() * nodes[0]->getNDofs();
	
    // temporary arrays and vector to hold the tangent, evectors and evalues
    matrix tg(n,n);
    matrix evectors(n,n);
    longvector evaluesR(n);
	longvector evaluesI(n);  // imaginary part.
	
    staticTangent(*this, tg);	
	tg.eigendata(evectors, evaluesR, evaluesI);
	
    Message("\n Element eigenvalues:");
    for (a=0; a<n; a++)
    {
        Message("\n     %2d : %+10.4e ", a+1, evaluesR.data[a]);
        im = evaluesI.data[a];
        if ( im*im  > 1e-12 )  Message(" %+10.4e *i", im );
    }
	
    Message("\n\n Element eigenvectors in the rows:");
	evectors.print();
}




void element ::	replaceNode(node& deleted, node& survivor)
{
	for (int i=0; i< getNNodes(); i++)
		if (nodes[i] == &deleted) nodes[i] = &survivor;
}




bool staticTangent(element& el, matrix& tangent)
{
	singleElementAssembler sea(el);
	sea._tasks.static_tangent = true;
    extern double ftgstiff;
    ftgstiff = 1.0;
	el.contributeToGenericIntegral(sea);
	tangent = sea.emat;
	return true;
}



bool element :: testImplementation() const
{
    return true;
}



