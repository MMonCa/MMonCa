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
 *  body.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 11/9/06.
 *  Modified by I.R. on 31/5/07
 *
 *
 * An important issue regarding the data in the body objects is the ownership of the 
 * element and node pointers. In the current design, the element and node pointers
 * in the object's data belong to this body. This means that, when the body is 
 * destructed, the data in these vectors must be deleted. When a nodeset or an elset
 * is created, the body lends the element/node pointers to this new class, but knowing
 * that if the nodeset/elset is destructed, the data associated with it should not
 * be destructed, because the body takes care of that.
 */

#include <algorithm>
#include <iomanip>
#include <iterator>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <map>
#include <set>
#include "boost/foreach.hpp"

#include "Model/Parts/body.h"
#include "General/feliksutil.h"

#include "Geometry/Solidmodel/BRep.h"


#include "Analysis/Integrators/integrator.h"
#include "Elements/evalspot.h"
#include "Io/logger.h"
#include "Io/message.h"
#include "Model/Sets/elset.h"
#include "Model/model.h"
#include "Model/Parts/modelpart.h"
#include "Math/tensor.h"
#include "Io/usercommand.h"

#ifdef WITHTBB
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/parallel_for_each.h"
#endif

using namespace blue;


body :: body(const string& name_, model& m, const size_t dim)
:	modelpart(name_, m),
	hasVariableMass(false),
	mass(0.0), massComputed(false),
	coordinationNumber(-33),
    _dimension(dim)
{	
}



body :: ~body()
{		
}


/* III

#ifdef WITHTBB
// data structure for tbb threaded version
class tbbAdvanceNodesInTime{
	public:
		tbbAdvanceNodesInTime( vector<node*>& nds, const integrator &integ, const double dt) 
			: my_nodes(&nds), my_integ(&integ), my_dt(dt){}

		void operator()( const tbb::blocked_range<size_t>& r) const
		{
			vector<node*> &nds( *my_nodes );
			for (size_t a=r.begin(); a!=r.end(); a++)
				nds[a]->advanceIntegration(*my_integ, my_dt);
			
		}
		
	private:
		vector<node*> *const my_nodes;
		integrator     const *my_integ;
		const double   my_dt;
};
#endif
 */



class advanceInTimeWorker
{
private:
    vector<node*> *const my_nodes;
    integrator     const *my_integ;
    const double   my_dt;
    
    
public:
    advanceInTimeWorker(vector<node*>& nds, const integrator &integ, const double dt) 
    : my_nodes(&nds), my_integ(&integ), my_dt(dt){}
    
    void operator()(node* nd) const { nd->advanceIntegration(*my_integ, my_dt);}
};




// when a solution has been obtained, all the body variables and the ones
// in its elements must be updated accordingly. 
void body :: advanceInTime(const integrator& integ, const double dt)
{
	if ( !isActive() ) return;

    updateCurrentState(dofset::tn1);
    feliks_for_each(elements.begin(),  elements.end(),  mem_fun(&element::commitCurrentState) );
    feliks_for_each(evalspots.begin(), evalspots.end(), mem_fun(&evalspot::commitCurrentState) );
    feliks_for_each(nodes.begin(),     nodes.end(),     advanceInTimeWorker(nodes, integ, dt));
    
	if (hasVariableMass)	computeMass();
}




bool body :: check()
{
	for_each(elements.begin(), elements.end(), mem_fun(&element::check) );
    for_each(evalspots.begin(), evalspots.end(), mem_fun(&evalspot::check) );
	return true;
}




void body :: computeExplicitAcceleration()
{
	feliks_for_each(nodes.begin(), nodes.end(), mem_fun(&node::computeExplicitAcceleration) );
}



void body :: computeExplicitEulerianAcceleration()
{
	feliks_for_each(nodes.begin(), nodes.end(), mem_fun(&node::computeExplicitEulerianAcceleration) );
}



void body :: computeExplicitVelocity()
{
	feliks_for_each(nodes.begin(), nodes.end(), mem_fun(&node::computeExplicitVelocity) );
}




const ivector body :: computeMaxSurfaceVelocity() const
{
	return getExternalNodes().getMaxVelocity();
}




void body :: computeMass()
{
	if ( !isActive() ) return;
    
    feliks_for_each(nodes.begin(),    nodes.end(),    mem_fun(&node::resetMass) );
    feliks_for_each(elements.begin(), elements.end(), mem_fun(&element::lumpMassToNodes) );
	feliks_for_each(evalspots.begin(), evalspots.end(), mem_fun(&evalspot::integrateNodalMasses) );
    
	mass = 0.0;
    BOOST_FOREACH(node* nd, nodes) mass += nd->getMass();

	massComputed = true;
    
}




void body :: computeReactionsSum(double *reac) const
{
	int dofs( linkedModel.getDofsPerNode() );
	
	Dzero(reac, dofs);
	vector<node*>::const_iterator iter = nodes.begin();
	node *ndp;
	while (iter != nodes.end() )
	{
		ndp = *iter;
		for (int i=0; i<dofs; i++) reac[i] += ndp->force()(i);
		++iter;
	}
}	




const ivector body :: getCurrentCenterOfMass()
{
	ivector com;
	if (mass == 0 || hasVariableMass) computeMass();
	
	// loop over all nodes, accumulating mass*position
	vector<node*>::iterator iter = nodes.begin();
	while ( iter != nodes.end() )
	{
		com += (*iter)->getMass() * (*iter)->coordinates();
		++iter;
	}
	
	com *= 1.0/mass;
	return com;
}



double body :: getMass() const
{
	if (!massComputed) const_cast<body*>(this)->computeMass();
	return mass;
}





#ifdef WITHTBB
// data structure for tbb threaded version
class tbbMaxElementEigenvalue{
public:
	double norm;
	
	tbbMaxElementEigenvalue( vector<element*>& elements) : my_elements(&elements), maxev(0.0), inElement(0){}
	tbbMaxElementEigenvalue( tbbMaxElementEigenvalue& nn, tbb::split) : my_elements(nn.my_elements), maxev(0.0), inElement(0){}
	void join(const tbbMaxElementEigenvalue& y)
    {
        if (maxev < y.maxev)
        {
            maxev = y.maxev;
            inElement = y.inElement;
        }
    }
	
	void operator()( const tbb::blocked_range<size_t>& r)
	{
		vector<element*> &nn( *my_elements );
		for (size_t a=r.begin(); a != r.end(); a++)
		{
			double elev = nn[a]->maxEigenvalue();
			if (elev > maxev)
			{
				maxev = elev;
				inElement = nn[a]->getLabel();
			}
		}
	}
	
public:
	vector<element*> *const my_elements;
	double	maxev;
	int		inElement;
};
#endif




/* finds the maximum eigenvalue for all the elements of the model. In fact, not quite
 the element eigenvalues as much as an upper bound.
 If an element returns 0 as maximum eigenvalue, it means that it does not compute any
 eigenvalue at all
 */
double body :: maxEigenvalueEstimate() const
{
    double  maxev(0.0);
	int maxLabel = 0;
	body* bb = const_cast<body*>(this);
	
#ifdef WITHTBB
    
	tbbMaxElementEigenvalue tbbmax( bb->elements );
	tbb::parallel_reduce(tbb::blocked_range<size_t>(0, elements.size()),
						 tbbmax,
						 tbb::auto_partitioner() );
	maxev    = tbbmax.maxev;
	maxLabel = tbbmax.inElement;
    
#else
	
	double eigv(0.0);
    BOOST_FOREACH(element* el, bb->elements)
	{
		eigv = el->maxEigenvalue();
		if (maxev < eigv)
		{
			maxev    = eigv;
			maxLabel = el->getLabel();
		}
	}
    
#endif
	
	if (maxev > 1e10)
	{
        std::cout << "enormous eigenvalue in element ";
		exit(0);		
	}
    return maxev;
}




feliks::meshing::solidmesh& body :: getTheSolidmesh()
{
    if ( theSolidmesh == 0 )
        throw runtime_error(" this body has no solidmesh associated to it");
    
   return *theSolidmesh;
}




const feliks::meshing::solidmesh& body :: getTheSolidmesh() const
{
    if ( theSolidmesh == 0 )
        throw runtime_error(" this body has no solidmesh associated to it");
    
    return *theSolidmesh;
}




set<node*> body :: getNodesFrom0Cells(const set<feliks::topology::cell*>& cells) const
{
    set<node*> s;
    
    BOOST_FOREACH(node* n, nodes)
    {
        feliks::topology::cell* c = const_cast<feliks::topology::cell*>(n->getCellP());
        if ( cells.find(c) != cells.end() ) 
            s.insert(n);
    }
    return s;
}




void body :: incrementSolution(const integrator& i)
{
	if (!active) return;
	
	if ( i.isNodallyBased() )
		i.incrementSolution( nodes );
	else
		i.incrementSolution( elements );    
}





void body :: initialize()
{
    if (isInitialized()) return;
    
    BOOST_FOREACH(element* el, elements)
    {
        el->setParentSubmesh(*this);
        el->initialize();
    }
    DebugMessage("Elements initialized");

    
    // Initialize nodes
	feliks_for_each( nodes.begin(), nodes.end(), mem_fun(&node::initializeDofs) );
	computeMass();
	DebugMessage("Nodes initialized");
		
	setInitialized();
	computeMass();
}




void body :: parse()
{
    
}




void body :: parseInformationReset()
{
    
}




void body :: print(std::ostream &os)
{
	os  << "\n Body name           : " << getName();
    modelpart::print(os);
}




void printBodyElements(const body& b, ostream &of)
{
	// check to avoid printing bodies with not elements
	if (b.elements.empty()) return;
	

	of << "\n\n  Body : " << b.getName();
	of << "\n  label  eltype   node1   node2   node3   node4 ...";

	vector<element*>::const_iterator iter=b.elements.begin();
	while (iter != b.elements.end())
	{
        of << endl << setw(6) << (*iter)->getLabel() << setw(6) << (*iter)->getEltype().label();
        for (int k=0; k< (*iter)->getNNodes(); k++)  
            of << setw(7) << (*iter)->getNode(k).getLabel();
		++iter;
	}
	of << "\n";
}



void printBodyFacesets(const body& b, ostream &of)
{
    /*
	if (facesets.size() == 0) return;
	
    std::map< std::string,faceset* >::const_iterator iter=b.getFacesets().begin();
	while ( iter != b.getFacesets().end() )
	{
		(*iter).second->info(of);
		++iter;
	}
	of << endl;
     */
}




void  printBodyReactionsInConstraints(const body& b)
{
	int dofsPerNode ( b.maxDofsPerNode() );
	double* ftot  = new double[dofsPerNode];
	Dzero(ftot, dofsPerNode);
	
    Message("\n\n           Nodal reactions (R = fint)");
    Message("\n          ---------------------------");
	
    vector<node*>::const_iterator iter = b.nodes.begin();
	node *ndp=NULL;
    while (iter != b.nodes.end() )
    {
		ndp = *iter;
        if (ndp->getNConstrainedDofs() > 0)
        {
            Message("\n Label(%4d) :", ndp->getLabel() );
            for (int i=0; i<dofsPerNode; i++)
            {
                if ( ndp->isDofConstrained(i) )
                {
                    Message("  % 11.5e", -ndp->force()(i) );
                    ftot[i] += ndp->force()(i);
                }
                else
                    Message("   -----------");
            }
        }
		++iter;
    }
    Message("\n          ----------------------------------");
    Message("\n Sum =        ");
    for (int i=0; i<dofsPerNode; i++) Message("  % 11.5e", -ftot[i]);
    delete []ftot;
}




void body :: remesh(bool& rec, bool& ghostnodes)
{
}



/*
void body :: replaceNode(node& deleted, node& survivor)
{	
	// in elements
	modelpart::replaceNode(deleted, survivor);
	
	// in nodesets
	getExternalNodes.replaceNode(deleted, survivor);
	allNodes.replaceNode(deleted, survivor);
	
	// delete the node in the global list
	nodes.erase( std::remove(nodes.begin(), nodes.end(), &deleted), nodes.end() );
	_nodesSorted = false;	
}
 */





// set all the nodal masses to zero
void  body :: resetNodalMasses()
{
	for_each( nodes.begin(), nodes.end(), mem_fun(&node::resetMass));	
}




// set all the nodal forces to zero
void  body :: resetNodalForces()
{
	feliks_for_each( nodes.begin(), nodes.end(), mem_fun(&node::forceReset));
}




// set all the nodal solution increments dU to zero
void  body :: resetSolutionIncrements()
{
	feliks_for_each( nodes.begin() , nodes.end(), mem_fun(&node::setDeltaDofsToZero) );
}




void body :: rewindSolution()
{
	// Recover nodal dofs from tn
	for_each( nodes.begin() , nodes.end(), mem_fun(&node::recoverFromBackupDofs) );
}




// this function scans a body from the commandline. Only regularbodies are defined in
// this way, but we put it here so that in the future maybe more elaborated things
// are added
body& body :: scan(const commandLine &cl, model &m)
{
	body *bd(0);	
	return *bd;
}




// one must be careful when renaming the body because the nodesets indexed with
// the body name are already in the map< >
void body :: setName(const string& newname)
{	
    modelpart::setName(newname);
	_externalNodes.setName ( newname + "_external_nodes");
	linkedModel.add(_externalNodes);
}




void body :: setNodalMasses(const longvector& v)
{
    vector<node*>::const_iterator iter = nodes.begin();
	node *ndp=NULL;
    while (iter != nodes.end() )
    {
		ndp = *iter;
		(*iter)->resetMass();
		(*iter)->incrementMass( v.data[ (*iter)->getID(0) ] );
		++iter;
    }	
}




void body :: updateCurrentState(const dofset::evaluation_time when)
{
    feliks_for_each(elements.begin(), elements.end(),
                  std::bind2nd(std::mem_fun(&element::updateCurrentState), when));
    
    feliks_for_each(evalspots.begin(), evalspots.end(),
                  std::bind2nd(std::mem_fun(&evalspot::updateCurrentState), when));   
}




void checkCoordNum(const double radioBusqueda, vector<body*> &bodies)
{
	double maxdiam(0.0);
    
	for (int i=1; i<bodies.size(); i++)
	{	
		bodies[i]->setNCoord(0);
		for (int j=i+1; j<bodies.size(); j++)
		{
			ivector cgi(bodies[i]->getCurrentCenterOfMass());
			ivector cgj(bodies[j]->getCurrentCenterOfMass());
			ivector rel(cgj-cgi);
			
			if (rel.norm() <= radioBusqueda)
			{
				bodies[i]->incrementNCoord();
				bodies[j]->incrementNCoord();
			}
		}
	}	
	
	extern double global_tn1;
	ofstream of("feliks.diameter", ofstream::app);
	of << setw(12) << scientific << global_tn1 << " " << maxdiam << "\n";
	of.close();
}

